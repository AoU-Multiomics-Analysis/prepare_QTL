"""Register and verify GitHub Apps without saving their credentials."""
import base64
import hmac
import html
import http.client
import http.server
import io
import json
import math
import re
import secrets
import time
import urllib.parse
import webbrowser

from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import padding, rsa

from release_setup import SetupError


REQUIRED_PERMISSIONS = {'contents': 'write', 'pull_requests': 'read'}


def _repository_owner(repository):
    if not isinstance(repository, str) or not re.fullmatch(
            r'[A-Za-z0-9][A-Za-z0-9-]{0,38}/[A-Za-z0-9_.-]{1,100}', repository):
        raise SetupError('Use a valid GitHub owner/repository name.')
    if repository.split('/')[1] in ('.', '..'):
        raise SetupError('Use a valid GitHub owner/repository name.')
    return repository.split('/')[0]


def app_manifest(repository: str, name: str, redirect_url: str) -> dict:
    _repository_owner(repository)
    if (not isinstance(name, str) or not 1 <= len(name) <= 100
            or any(ord(c) < 32 for c in name)):
        raise SetupError('Use a valid GitHub App name.')
    if not isinstance(redirect_url, str) or not re.fullmatch(
            r'http://127\.0\.0\.1:[0-9]{1,5}/callback', redirect_url):
        raise SetupError('The App callback must use the local loopback address.')
    port = int(redirect_url.split(':')[2].split('/')[0])
    if not 1 <= port <= 65535:
        raise SetupError('The App callback must use the local loopback address.')
    return {
        'name': name, 'url': f'https://github.com/{repository}',
        'redirect_url': redirect_url, 'public': False,
        'hook_attributes': {'url': f'https://github.com/{repository}', 'active': False},
        'default_permissions': dict(REQUIRED_PERMISSIONS), 'default_events': [],
    }


def validate_callback(path: str, query: str, expected_state: str) -> str:
    error = 'The App callback was rejected.'
    if (path != '/callback' or not isinstance(query, str) or len(query) > 2048
            or not isinstance(expected_state, str) or not expected_state):
        raise SetupError(error)
    try:
        pairs = urllib.parse.parse_qsl(query, keep_blank_values=True,
                                      strict_parsing=True, max_num_fields=2)
        fields = dict(pairs)
        valid_state = hmac.compare_digest(fields.get('state', ''), expected_state)
    except (ValueError, TypeError):
        raise SetupError(error) from None
    if (len(pairs) != 2 or set(fields) != {'state', 'code'} or not valid_state
            or not re.fullmatch(r'[A-Za-z0-9_-]{1,512}', fields['code'])):
        raise SetupError(error)
    return fields['code']


def app_jwt(app_id: int, pem: bytes, now: int) -> str:
    if type(app_id) is not int or app_id <= 0 or type(now) is not int:
        raise SetupError('Use a valid GitHub App ID and RSA private key.')
    try:
        key = serialization.load_pem_private_key(pem, password=None)
        if not isinstance(key, rsa.RSAPrivateKey) or key.key_size < 2048:
            raise ValueError()

        def encode(value):
            return base64.urlsafe_b64encode(value).rstrip(b'=')

        header = encode(b'{"alg":"RS256","typ":"JWT"}')
        payload = encode(json.dumps({'iat': now - 60, 'exp': now + 540,
                                     'iss': str(app_id)}, separators=(',', ':')).encode())
        message = header + b'.' + payload
        signature = key.sign(message, padding.PKCS1v15(), hashes.SHA256())
        return (message + b'.' + encode(signature)).decode('ascii')
    except Exception:
        raise SetupError('Use a valid GitHub App ID and RSA private key.') from None


def _permissions(value):
    if not isinstance(value, dict) or any(
            not isinstance(name, str) or not re.fullmatch(r'[a-z_]{1,100}', name)
            or access not in ('read', 'write', 'admin')
            for name, access in value.items()):
        raise SetupError('GitHub returned invalid App permissions.')
    if (value.get('contents') != 'write'
            or value.get('pull_requests') not in ('read', 'write')):
        raise SetupError('The App needs Contents write and Pull requests read permissions.')
    return dict(value)


def _app_identity(value, app_id, owner):
    if (not isinstance(value, dict) or type(value.get('id')) is not int
            or value['id'] != app_id
            or not isinstance(value.get('slug'), str)
            or not re.fullmatch(r'[a-z0-9][a-z0-9-]{0,99}', value['slug'])
            or not isinstance(value.get('owner'), dict)
            or str(value['owner'].get('login', '')).lower() != owner.lower()):
        raise SetupError('The GitHub App identity or owner does not match.')
    return value['slug']


def verify_installation(repository: str, app_id: int, pem: bytes, *,
                        newly_created: bool = False) -> dict:
    """Read App access without creating an installation token or changing scope.

    The caller sets newly_created only for this run's manifest registration.
    Notices and all returned fields are public; the JWT is never returned.
    Selected mode proves target access, not the absence of other selected repos.
    """
    owner = _repository_owner(repository)
    token = app_jwt(app_id, pem, int(time.time()))
    try:
        identity = _github_json('GET', '/app', token=token)
        slug = _app_identity(identity, app_id, owner)
        requested_permissions = _permissions(identity.get('permissions'))
        installed = _github_json('GET', f'/repos/{repository}/installation', token=token)
        if (not isinstance(installed, dict) or type(installed.get('app_id')) is not int
                or installed['app_id'] != app_id or installed.get('app_slug') != slug
                or not isinstance(installed.get('account'), dict)
                or str(installed['account'].get('login', '')).lower() != owner.lower()):
            raise SetupError('The App installation identity or owner does not match.')
        if (installed.get('suspended_at', True) is not None
                or installed.get('suspended_by') is not None):
            raise SetupError('The App installation is suspended or has an invalid status.')
        permissions = _permissions(installed.get('permissions'))
        extra_write = any(access != 'read' and (name, access) != ('contents', 'write')
                          for entry in (requested_permissions, permissions)
                          for name, access in entry.items())
        selection = installed.get('repository_selection')
        if selection not in ('selected', 'all'):
            raise SetupError('GitHub returned an invalid repository selection.')
        if newly_created and (extra_write or selection == 'all'):
            raise SetupError('The new App has broader access than required; review its settings.')
        notices = []
        if selection == 'all':
            notices.append('The existing App installation can access all owner repositories.')
        else:
            notices.append('Target repository access is verified; other selected repositories were not checked.')
        if extra_write:
            notices.append('The existing App has additional write permissions; access was not changed.')
        return {'id': app_id, 'slug': slug, 'repository': repository,
                'permissions': permissions, 'notices': notices}
    except SetupError:
        raise
    except Exception:
        raise SetupError('GitHub App verification failed; no response details were retained.') from None


def _github_json(method, path, *, token=None):
    """Use direct TLS to GitHub; redirects, proxies and token creation are absent."""
    allowed = (method == 'GET' and (path == '/app' or re.fullmatch(
        r'/repos/[A-Za-z0-9][A-Za-z0-9-]{0,38}/[A-Za-z0-9_.-]{1,100}/installation', path)))
    allowed = allowed or (method == 'POST' and re.fullmatch(
        r'/app-manifests/[A-Za-z0-9_-]{1,512}/conversions', path))
    if not allowed:
        raise SetupError('The GitHub App API request was rejected.')
    connection = None
    try:
        connection = http.client.HTTPSConnection('api.github.com', timeout=15)
        headers = {'Accept': 'application/vnd.github+json',
                   'X-GitHub-Api-Version': '2022-11-28', 'User-Agent': 'release-setup'}
        if token is not None:
            headers['Authorization'] = f'Bearer {token}'
        connection.request(method, path, headers=headers)
        response = connection.getresponse()
        if not 200 <= response.status < 300:
            raise ValueError()
        payload = response.read(1048577)
        if len(payload) > 1048576:
            raise ValueError()
        result = json.loads(payload)
        if not isinstance(result, dict):
            raise ValueError()
        return result
    except Exception:
        raise SetupError('GitHub App API access failed; no response details were retained.') from None
    finally:
        if connection is not None:
            try:
                connection.close()
            except Exception:
                raise SetupError('The GitHub App API connection could not close.') from None


class _DeadlineReader(io.RawIOBase):
    """Limit a request by elapsed time and bytes, including slowly sent headers."""

    def __init__(self, connection, deadline):
        self.connection = connection
        self.deadline = deadline
        self.remaining = 16384

    def readable(self):
        return True

    def readinto(self, buffer):
        remaining_time = self.deadline - time.monotonic()
        if remaining_time <= 0 or self.remaining <= 0:
            raise TimeoutError()
        self.connection.settimeout(remaining_time)
        size = self.connection.recv_into(buffer, min(len(buffer), self.remaining))
        self.remaining -= size
        return size


class _CallbackServer(http.server.HTTPServer):
    def handle_error(self, request, client_address):
        # HTTPServer normally prints exception details, which can contain a code.
        pass


def _callback_code(repository, owner_type, timeout):
    state = secrets.token_urlsafe(32)
    start_path = '/start/' + secrets.token_urlsafe(32)
    deadline = time.monotonic() + timeout
    code = None

    class Handler(http.server.BaseHTTPRequestHandler):
        def setup(self):
            super().setup()
            self.rfile.close()
            self.rfile = io.BufferedReader(_DeadlineReader(
                self.connection, min(deadline, time.monotonic() + 2)))

        def log_message(self, format, *args):
            pass

        def send_error(self, code, message=None, explain=None):
            self.respond(400, b'The request was rejected.')

        def respond(self, status, body):
            self.close_connection = True
            self.send_response(status)
            self.send_header('Content-Type', 'text/html; charset=utf-8')
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('Referrer-Policy', 'no-referrer')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Content-Security-Policy',
                             "default-src 'none'; form-action https://github.com; "
                             "frame-ancestors 'none'; base-uri 'none'")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            nonlocal code
            if (self.headers.get_all('Host') != [host]
                    or self.path != self.requestline.split()[1]
                    or time.monotonic() >= deadline):
                self.respond(400, b'The request was rejected.')
                return
            if self.path == start_path:
                self.respond(200, page)
                return
            # Split only the exact origin-form target; reject absolute URLs.
            path, _, query = self.path.partition('?')
            if path != '/callback':
                self.respond(404, b'The page was not found.')
                return
            try:
                validated = validate_callback(path, query, state)
            except SetupError:
                self.respond(400, b'The callback was rejected.')
                return
            code = validated
            self.respond(200, b'App registration received. Return to the terminal.')

    with _CallbackServer(('127.0.0.1', 0), Handler) as server:
        host = f'127.0.0.1:{server.server_port}'
        callback = f'http://{host}/callback'
        manifest = app_manifest(repository, 'Release helper', callback)
        owner = _repository_owner(repository)
        route = (f'/organizations/{owner}/settings/apps/new'
                 if owner_type == 'Organization' else '/settings/apps/new')
        action = 'https://github.com' + route + '?' + urllib.parse.urlencode({'state': state})
        page = ('<!doctype html><html><head><meta charset="utf-8">'
                '<title>Register release App</title></head><body>'
                '<p>Continue to GitHub to register the release App.</p>'
                '<form method="post" action="' + html.escape(action, quote=True) + '">'
                '<input type="hidden" name="manifest" value="'
                + html.escape(json.dumps(manifest, separators=(',', ':')), quote=True)
                + '"><button type="submit">Continue to GitHub</button></form></body></html>').encode()
        if not webbrowser.open(f'http://{host}{start_path}'):
            raise SetupError('The browser could not open the App registration page.')
        while code is None and time.monotonic() < deadline:
            server.timeout = min(0.2, max(0, deadline - time.monotonic()))
            server.handle_request()
        if code is None:
            raise SetupError('App registration timed out; check GitHub before trying again.')
    return code


def register_app(repository: str, owner_type: str, *, timeout: int = 600) -> dict:
    """Open the consent page and return credentials only in memory.

    The CLI displays the public ID and installation URL, waits for the user to
    install on the target repository, then calls verify_installation.
    """
    owner = _repository_owner(repository)
    if owner_type not in ('Organization', 'User'):
        raise SetupError('GitHub App registration needs a user or organization owner.')
    if (not isinstance(timeout, (int, float)) or isinstance(timeout, bool)
            or not math.isfinite(timeout) or not 0 < timeout <= 600):
        raise SetupError('Use an App registration timeout between zero and 600 seconds.')
    try:
        code = _callback_code(repository, owner_type, timeout)
        converted = _github_json('POST', f'/app-manifests/{code}/conversions')
        app_id = converted.get('id')
        if type(app_id) is not int or app_id <= 0:
            raise ValueError()
        slug = _app_identity(converted, app_id, owner)
        permissions = _permissions(converted.get('permissions'))
        if any(access != 'read' and (name, access) != ('contents', 'write')
               for name, access in permissions.items()):
            raise SetupError('The new App has broader access than required; review its settings.')
        pem = converted['pem'].encode('ascii')
        app_jwt(app_id, pem, int(time.time()))
        return {'id': app_id, 'slug': slug, 'pem': pem}
    except SetupError:
        raise
    except Exception:
        raise SetupError('App registration failed; check GitHub before trying again.') from None
