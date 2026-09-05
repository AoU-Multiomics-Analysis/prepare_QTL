"""Guard the local App callback and in-memory credentials."""
import base64
import contextlib
import io
import http.client
import json
import os
import socket
import sys
import threading
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path
import unittest
from unittest import mock

from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import padding, rsa

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ci'))
from release_setup import SetupError
import release_app_setup as app

SECRET = 'FAKE-SECRET-MUST-NOT-ESCAPE'


class CLIKeyFileTests(unittest.TestCase):
    def test_private_regular_file_is_read_and_unsafe_inputs_are_redacted(self):
        from setup_release import _read_private_key
        with tempfile.TemporaryDirectory() as directory:
            key_path = Path(directory) / 'test.pem'
            key_path.write_bytes(SECRET.encode())
            key_path.chmod(0o600)
            self.assertEqual(_read_private_key(str(key_path)), SECRET.encode())
            for path in (directory, str(key_path) + '.missing'):
                with self.subTest(path=path), self.assertRaises(SetupError) as caught:
                    _read_private_key(path)
                self.assertNotIn(SECRET, str(caught.exception))
            if os.name == 'posix':
                key_path.chmod(0o640)
                with self.assertRaises(SetupError) as caught:
                    _read_private_key(str(key_path))
                self.assertIn('chmod 600', str(caught.exception))
                self.assertEqual(key_path.stat().st_mode & 0o777, 0o640)
                key_path.chmod(0o604)
                with self.assertRaises(SetupError):
                    _read_private_key(str(key_path))

    @unittest.skipUnless(os.name == 'posix', 'POSIX special files')
    def test_symlink_and_fifo_are_rejected_without_blocking(self):
        from setup_release import _read_private_key
        with tempfile.TemporaryDirectory() as directory:
            key_path = Path(directory) / 'test.pem'
            key_path.write_bytes(SECRET.encode())
            key_path.chmod(0o600)
            link = Path(directory) / 'link.pem'
            link.symlink_to(key_path)
            fifo = Path(directory) / 'fifo.pem'
            os.mkfifo(fifo, 0o600)
            for path in (link, fifo):
                with self.subTest(path=path), self.assertRaises(SetupError):
                    _read_private_key(str(path))


def public_app():
    return {'id': 42, 'slug': 'release-helper', 'owner': {'login': 'owner'},
            'permissions': {'contents': 'write', 'pull_requests': 'read',
                            'metadata': 'read'}}


def installation():
    return {'id': 7, 'app_id': 42, 'app_slug': 'release-helper',
            'account': {'login': 'owner'}, 'suspended_at': None,
            'suspended_by': None, 'repository_selection': 'selected',
            'permissions': {'contents': 'write', 'pull_requests': 'read',
                            'metadata': 'read'}}


class CredentialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
        cls.pem = cls.key.private_bytes(serialization.Encoding.PEM,
                                       serialization.PrivateFormat.PKCS8,
                                       serialization.NoEncryption())

    def test_jwt_claims_expiry_and_signature(self):
        token = app.app_jwt(42, self.pem, 1700000000)
        header, payload, signature = token.split('.')
        decode = lambda value: base64.urlsafe_b64decode(value + '=' * (-len(value) % 4))
        self.assertEqual(json.loads(decode(header)), {'alg': 'RS256', 'typ': 'JWT'})
        self.assertEqual(json.loads(decode(payload)),
                         {'iat': 1699999940, 'exp': 1700000540, 'iss': '42'})
        self.assertNotIn('=', token)
        self.key.public_key().verify(decode(signature), f'{header}.{payload}'.encode(),
                                     padding.PKCS1v15(), hashes.SHA256())

    def test_invalid_credentials_are_redacted(self):
        for app_id, pem in [(42, SECRET.encode()), (True, self.pem), (0, self.pem)]:
            with self.assertRaises(SetupError) as error:
                app.app_jwt(app_id, pem, 1700000000)
            self.assertNotIn(SECRET, str(error.exception))
            self.assertTrue(error.exception.__suppress_context__ or pem == self.pem)

    def verify(self, first=None, second=None, *, newly_created=False):
        responses = [first if first is not None else public_app(),
                     second if second is not None else installation()]
        calls = []

        def api(method, path, *, token=None):
            calls.append((method, path, token))
            return responses.pop(0)

        with mock.patch.object(app, '_github_json', api, create=True):
            result = app.verify_installation('owner/repo', 42, self.pem,
                                             newly_created=newly_created)
        self.assertEqual([(m, p) for m, p, _ in calls],
                         [('GET', '/app'), ('GET', '/repos/owner/repo/installation')])
        self.assertTrue(all(token.count('.') == 2 for _, _, token in calls))
        return result

    def test_verification_returns_public_fields_only_and_never_mints_token(self):
        first, second = public_app(), installation()
        first['client_secret'] = second['sensitive'] = SECRET
        result = self.verify(first, second)
        self.assertEqual(result['id'], 42)
        self.assertEqual(result['slug'], 'release-helper')
        self.assertEqual(result['repository'], 'owner/repo')
        self.assertEqual(result['permissions'],
                         {'contents': 'write', 'pull_requests': 'read', 'metadata': 'read'})
        self.assertEqual(set(result), {'id', 'slug', 'repository', 'permissions', 'notices'})
        self.assertNotIn(SECRET, json.dumps(result))

    def test_verification_rejects_identity_owner_permissions_and_suspension(self):
        cases = [
            ('app', 'id', 99), ('app', 'owner', {'login': 'other'}),
            ('app', 'slug', SECRET + '\n'),
            ('install', 'app_id', 99), ('install', 'app_slug', 'other'),
            ('install', 'account', {'login': 'other'}),
            ('install', 'suspended_at', '2026-09-05T00:00:00Z'),
            ('install', 'suspended_by', {'login': 'admin'}),
            ('install', 'repository_selection', 'unknown'),
        ]
        for which, field, value in cases:
            first, second = public_app(), installation()
            (first if which == 'app' else second)[field] = value
            with self.subTest(which=which, field=field), self.assertRaises(SetupError):
                self.verify(first, second)
        for which in ('app', 'install'):
            for permissions in ({'contents': 'read', 'pull_requests': 'read'},
                                {'contents': 'write'}):
                first, second = public_app(), installation()
                (first if which == 'app' else second)['permissions'] = permissions
                with self.subTest(which=which, permissions=permissions), self.assertRaises(SetupError):
                    self.verify(first, second)

    def test_broad_existing_installation_is_reported_but_new_app_is_rejected(self):
        second = installation()
        second['repository_selection'] = 'all'
        self.assertTrue(self.verify(second=second)['notices'])
        with self.assertRaises(SetupError):
            self.verify(second=second, newly_created=True)
        first = public_app()
        first['permissions']['issues'] = 'write'
        self.assertTrue(self.verify(first=first)['notices'])
        with self.assertRaises(SetupError):
            self.verify(first=first, newly_created=True)

    def test_transport_failure_cannot_leak_error_body(self):
        output, errors = io.StringIO(), io.StringIO()
        with mock.patch.object(app, '_github_json', side_effect=RuntimeError(SECRET), create=True):
            with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
                with self.assertRaises(SetupError) as raised:
                    app.verify_installation('owner/repo', 42, self.pem)
        self.assertNotIn(SECRET, output.getvalue() + errors.getvalue() + str(raised.exception))


class CallbackTests(unittest.TestCase):
    def test_callback_requires_exact_state_path_and_single_code(self):
        self.assertEqual(app.validate_callback('/callback', 'state=s&code=c', 's'), 'c')
        for path, query in [
            ('/other', 'state=s&code=c'), ('/callback', 'state=wrong&code=c'),
            ('/callback', 'state=s&code=c&code=d'), ('/callback', 'state=s'),
            ('/callback', 'state=s&state=s&code=c'),
            ('/callback', 'state=s&code=c&extra=x'),
            ('/callback', 'state=s&code=../bad'),
            ('/callback', 'state=s&code=%0Asecret'),
            ('/callback', 'state=s&code='),
            ('/callback', 'state=s&code=' + 'a' * 513),
            ('/callback', 'state=s&code=%ZZ'),
        ]:
            with self.subTest(path=path, query=query), self.assertRaises(SetupError):
                app.validate_callback(path, query, 's')

    def test_manifest_has_only_required_permissions(self):
        manifest = app.app_manifest('owner/repo', 'Release "helper"',
                                    'http://127.0.0.1:12345/callback')
        self.assertEqual(manifest['default_permissions'],
                         {'contents': 'write', 'pull_requests': 'read'})
        self.assertEqual(manifest['default_events'], [])
        self.assertFalse(manifest['hook_attributes']['active'])
        self.assertFalse(manifest['public'])
        self.assertEqual(manifest['url'], 'https://github.com/owner/repo')
        self.assertEqual(manifest['redirect_url'], 'http://127.0.0.1:12345/callback')

    def test_manifest_rejects_nonlocal_redirect_and_unsafe_repository(self):
        for repo, url in [('owner/repo', 'https://evil.test/callback'),
                          ('owner/repo', 'http://127.0.0.1:1/other'),
                          ('owner/repo', 'http://127.0.0.1:99999/callback'),
                          ('owner/repo/../x', 'http://127.0.0.1:12345/callback')]:
            with self.assertRaises(SetupError):
                app.app_manifest(repo, 'helper', url)


class TransportTests(unittest.TestCase):
    def test_close_failure_cannot_expose_secret(self):
        connection = mock.MagicMock()
        connection.getresponse.side_effect = TimeoutError(SECRET)
        connection.close.side_effect = OSError(SECRET)
        with mock.patch('http.client.HTTPSConnection', return_value=connection):
            with self.assertRaises(SetupError) as raised:
                app._github_json('GET', '/app', token=SECRET)
        self.assertNotIn(SECRET, str(raised.exception))

    def test_transport_uses_fixed_https_host_and_does_not_follow_redirects(self):
        observations = []

        class Connection:
            status = 200

            def __init__(self, host, *, timeout):
                observations.append(('connect', host, timeout))

            def request(self, method, path, *, body=None, headers):
                observations.append(('request', method, path, headers))

            def getresponse(self):
                return self

            def read(self, limit):
                return b'{"id":42}'

            def close(self):
                observations.append(('close',))

        with mock.patch('http.client.HTTPSConnection', Connection):
            self.assertEqual(app._github_json('GET', '/app', token='jwt'), {'id': 42})
            self.assertEqual(observations[0][1], 'api.github.com')
            self.assertGreater(observations[0][2], 0)
            self.assertEqual(observations[1][1:3], ('GET', '/app'))
            self.assertEqual(observations[1][3]['Authorization'], 'Bearer jwt')
            self.assertEqual(observations[-1], ('close',))
            for status in (301, 302, 307, 308, 401, 403, 500):
                Connection.status = status
                observations.clear()
                with self.assertRaises(SetupError):
                    app._github_json('POST', '/app-manifests/temporary/conversions')
                self.assertEqual(len([x for x in observations if x[0] == 'connect']), 1)
                self.assertEqual(observations[-1], ('close',))

    def test_transport_redacts_network_json_and_oversized_response_errors(self):
        for outcome in [TimeoutError(SECRET), ValueError(SECRET),
                        SECRET.encode(), b'[' + b'0,' * 600000 + b'0]']:
            connection = mock.MagicMock()
            connection.getresponse.return_value.status = 200
            if isinstance(outcome, Exception):
                connection.getresponse.side_effect = outcome
            else:
                connection.getresponse.return_value.read.return_value = outcome
            with mock.patch('http.client.HTTPSConnection', return_value=connection):
                with self.assertRaises(SetupError) as raised:
                    app._github_json('GET', '/app', token=SECRET)
            self.assertNotIn(SECRET, str(raised.exception))
            self.assertTrue(raised.exception.__suppress_context__)

    def test_transport_rejects_unapproved_paths(self):
        for method, path in [('GET', 'https://evil.test/app'),
                             ('POST', '/app/installations/7/access_tokens'),
                             ('GET', '/repos/owner/repo/installation?evil=yes'),
                             ('POST', '/app-manifests/../../evil/conversions')]:
            with self.assertRaises(SetupError):
                app._github_json(method, path, token=SECRET)


class RegistrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
        cls.pem = key.private_bytes(serialization.Encoding.PEM,
                                   serialization.PrivateFormat.PKCS8,
                                   serialization.NoEncryption())

    def conversion(self):
        return {**public_app(), 'pem': self.pem.decode(),
                'client_secret': SECRET, 'webhook_secret': SECRET}

    def run_browser(self, interaction, *, owner_type='Organization', conversion=None,
                    timeout=3):
        workers, failures, urls = [], [], []

        def open_browser(url):
            urls.append(url)

            def run():
                try:
                    interaction(url)
                except Exception as error:
                    failures.append(error)

            worker = threading.Thread(target=run, daemon=True)
            workers.append(worker)
            worker.start()
            return True

        converted = []

        def exchange(method, path, *, token=None):
            converted.append((method, path, token))
            if isinstance(conversion, Exception):
                raise conversion
            return self.conversion() if conversion is None else conversion

        output, errors = io.StringIO(), io.StringIO()
        try:
            with mock.patch('webbrowser.open', open_browser), \
                    mock.patch.object(app, '_github_json', exchange), \
                    contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
                result = app.register_app('owner/repo', owner_type, timeout=timeout)
        finally:
            for worker in workers:
                worker.join(4)
                self.assertFalse(worker.is_alive())
            self.assertNotIn(SECRET, output.getvalue() + errors.getvalue())
            if failures:
                raise failures[0]
            if urls:
                parsed = urllib.parse.urlsplit(urls[0])
                with socket.socket() as client:
                    client.settimeout(0.2)
                    self.assertNotEqual(client.connect_ex(('127.0.0.1', parsed.port)), 0)
        return result, converted

    def fetch_form(self, url, owner_type='Organization'):
        with urllib.request.urlopen(url, timeout=2) as response:
            html = response.read().decode()
            self.assertEqual(response.headers['Cache-Control'], 'no-store')
            self.assertEqual(response.headers['Referrer-Policy'], 'no-referrer')
            self.assertIn("default-src 'none'", response.headers['Content-Security-Policy'])
        from html.parser import HTMLParser

        class FormParser(HTMLParser):
            def handle_starttag(self, tag, attributes):
                fields = dict(attributes)
                if tag == 'form':
                    self.form = fields
                elif tag == 'input' and fields.get('name') == 'manifest':
                    self.manifest = json.loads(fields['value'])

        parser = FormParser()
        parser.feed(html)
        self.assertEqual(parser.form['method'].lower(), 'post')
        action = urllib.parse.urlsplit(parser.form['action'])
        self.assertEqual(action.netloc, 'github.com')
        self.assertEqual(action.scheme, 'https')
        self.assertEqual(action.path, '/organizations/owner/settings/apps/new'
                         if owner_type == 'Organization' else '/settings/apps/new')
        self.assertEqual(parser.manifest['default_permissions'],
                         {'contents': 'write', 'pull_requests': 'read'})
        state = urllib.parse.parse_qs(action.query)['state'][0]
        self.assertGreaterEqual(len(state), 40)
        self.assertNotEqual(urllib.parse.urlsplit(url).path, '/callback')
        return parser.manifest['redirect_url'], state

    def finish(self, callback, state):
        query = urllib.parse.urlencode({'state': state, 'code': 'temporary-code'})
        with urllib.request.urlopen(callback + '?' + query, timeout=2) as response:
            self.assertEqual(response.status, 200)

    def test_org_and_personal_flow_closes_server_and_discards_other_secrets(self):
        for owner_type in ('Organization', 'User'):
            def browser(url):
                callback, state = self.fetch_form(url, owner_type)
                self.finish(callback, state)

            result, calls = self.run_browser(browser, owner_type=owner_type)
            self.assertEqual(result, {'id': 42, 'slug': 'release-helper', 'pem': self.pem})
            self.assertEqual(calls, [('POST', '/app-manifests/temporary-code/conversions', None)])

    def test_host_unknown_path_and_bad_callback_do_not_complete_registration(self):
        def browser(url):
            parsed = urllib.parse.urlsplit(url)
            for path, host in [(parsed.path, 'evil.test'), ('/unknown', parsed.netloc),
                               ('/callback?state=wrong&code=' + SECRET, parsed.netloc)]:
                client = http.client.HTTPConnection('127.0.0.1', parsed.port, timeout=2)
                client.request('GET', path, headers={'Host': host})
                response = client.getresponse()
                body = response.read().decode()
                self.assertIn(response.status, (400, 404))
                self.assertNotIn(SECRET, body)
                client.close()
            callback, state = self.fetch_form(url)
            self.finish(callback, state)

        result, calls = self.run_browser(browser)
        self.assertEqual(result['id'], 42)
        self.assertEqual(len(calls), 1)

    def test_timeout_including_partial_request_closes_server(self):
        started = time.monotonic()

        def browser(url):
            parsed = urllib.parse.urlsplit(url)
            with socket.create_connection(('127.0.0.1', parsed.port), timeout=2) as client:
                client.sendall(b'GET /callback?code=' + SECRET.encode())
                while client.recv(4096):
                    pass

        with self.assertRaises(SetupError):
            self.run_browser(browser, timeout=0.3)
        self.assertLess(time.monotonic() - started, 2)

    def test_conversion_failure_has_fixed_error_without_sensitive_output(self):
        def browser(url):
            self.finish(*self.fetch_form(url))

        for conversion in [RuntimeError(SECRET), {'pem': SECRET},
                           {**self.conversion(), 'owner': {'login': 'other'}}]:
            with self.assertRaises(SetupError) as raised:
                self.run_browser(browser, conversion=conversion)
            self.assertNotIn(SECRET, str(raised.exception))

    def test_double_slash_callback_is_rejected_before_valid_callback(self):
        def browser(url):
            callback, state = self.fetch_form(url)
            parsed = urllib.parse.urlsplit(callback)
            client = http.client.HTTPConnection('127.0.0.1', parsed.port, timeout=2)
            client.request('GET', '//callback?' + urllib.parse.urlencode(
                {'state': state, 'code': 'wrong-code'}))
            response = client.getresponse()
            self.assertEqual(response.status, 400)
            response.read()
            client.close()
            self.finish(callback, state)

        _, calls = self.run_browser(browser)
        self.assertEqual(calls, [('POST', '/app-manifests/temporary-code/conversions', None)])

    def test_slow_header_bytes_cannot_extend_request_deadline(self):
        started = time.monotonic()

        def browser(url):
            parsed = urllib.parse.urlsplit(url)
            with socket.create_connection(('127.0.0.1', parsed.port), timeout=1) as client:
                client.sendall(b'GET / HTTP/1.1\r\nHost: ')
                while time.monotonic() - started < 1.5:
                    try:
                        client.sendall(b'a')
                    except OSError:
                        break
                    time.sleep(0.03)

        with self.assertRaises(SetupError):
            self.run_browser(browser, timeout=0.3)
        self.assertLess(time.monotonic() - started, 1)
