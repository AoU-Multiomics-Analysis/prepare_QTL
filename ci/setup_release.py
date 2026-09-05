#!/usr/bin/env python3
"""Inspect release settings, or explicitly apply missing settings and credentials."""
import argparse
import os
import re
import stat

from release_setup import (GitHubCLI, SetupError, apply_settings,
                           inspect_settings, plan_settings)


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        # argparse normally repeats unknown argument values, which could be keys.
        super().error('Invalid command arguments. Use --help; private keys require --private-key-file.')


def _arguments(argv):
    parser = _Parser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--repo', required=True, metavar='OWNER/REPO')
    parser.add_argument('--reviewer', required=True, action='append', metavar='LOGIN',
                        help='Required reviewer; repeat for one to six unique GitHub users.')
    parser.add_argument('--apply', action='store_true', help='Create missing settings. Default: read-only plan.')
    parser.add_argument('--create-app', action='store_true', help='Guide new App registration and installation.')
    parser.add_argument('--app-id', type=int, help='Existing GitHub App ID; requires a private-key file.')
    parser.add_argument('--private-key-file', metavar='PATH', help='Private regular file; mode 600 on POSIX. Never read in dry run.')
    args = parser.parse_args(argv)
    if (not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9-]{0,38}/[A-Za-z0-9_.-]{1,100}', args.repo)
            or args.repo.split('/')[1] in ('.', '..')):
        parser.error('Invalid repository.')
    if (not 1 <= len(args.reviewer) <= 6
            or len({name.lower() for name in args.reviewer}) != len(args.reviewer)
            or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9-]{0,38}', name)
                   for name in args.reviewer)):
        parser.error('Invalid reviewers.')
    if ((args.app_id is None) != (args.private_key_file is None)
            or args.app_id is not None and args.app_id <= 0
            or args.create_app and args.app_id is not None):
        parser.error('Invalid credential mode.')
    return args


def _read_private_key(path):
    """Check the opened file before reading, without changing its permissions."""
    descriptor = None
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0)
                             | getattr(os, 'O_NONBLOCK', 0))
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise SetupError('The private-key path must be a regular file.')
        if os.name == 'posix' and metadata.st_mode & 0o077:
            raise SetupError('The private-key file is accessible to other users. Run chmod 600 on the file, then retry.')
        with os.fdopen(descriptor, 'rb') as stream:
            descriptor = None
            return stream.read()
    except OSError:
        raise SetupError('The private-key file could not be opened or read. Use a regular file, not a symbolic link.') from None
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _unchanged(expected, current):
    # Registry availability is a notice, not a mutable repository setting.
    if ({key: value for key, value in expected.items() if key != 'readiness_notices'}
            != {key: value for key, value in current.items() if key != 'readiness_notices'}):
        raise SetupError('Settings changed after inspection; stop other setup operators and inspect again.')


def _verify_settings(snapshot, app_id):
    plan = plan_settings(snapshot, app_id)
    if plan['conflicts']:
        for conflict in plan['conflicts']:
            print(f'Conflict: {conflict}')
        raise SetupError('Release settings could not be verified. Resolve the listed conflicts.')
    if plan['actions']:
        raise SetupError('Required release settings are missing; inspect and rerun setup.')


def main(argv=None, *, gh=None, app_service=None) -> int:
    args = _arguments(argv)
    gh = GitHubCLI() if gh is None else gh
    stage = 'inspect'
    completed = []
    app_id = args.app_id
    pem = None
    uploaded = False
    upload_attempted = False
    try:
        print('stage=inspect status=started')
        snapshot = inspect_settings(gh, args.repo, args.reviewer)
        plan = plan_settings(snapshot, args.app_id)
        if args.create_app and (snapshot['secret_present']
                                or 'RELEASE_APP_ID' in snapshot['variables']):
            plan['conflicts'].append(
                'App credentials already exist. Use --app-id and --private-key-file for the existing App.')
        if args.apply and not args.create_app and args.app_id is None:
            plan['conflicts'].append(
                'Apply requires --create-app or both --app-id and --private-key-file for installation verification.')
        print('stage=inspect status=complete')
        stage = 'plan'
        if not plan['conflicts']:
            for action in plan['actions']:
                print(f"Create {action['kind']}:{action['target']}")
            if args.create_app:
                print('Register an App, then wait for installation approval.')
            elif args.app_id is not None:
                print(f'Verify supplied credentials for App ID: {args.app_id}')
            if not snapshot['secret_present']:
                print('Upload an absent RELEASE_APP_PRIVATE_KEY only after App and environment verification.')
            if 'RELEASE_APP_ID' not in snapshot['variables']:
                print('Create absent RELEASE_APP_ID only after key upload.')
        for notice in plan['notices']:
            print(f'Notice: {notice}')
        for conflict in plan['conflicts']:
            print(f'Conflict: {conflict}')
        if plan['conflicts']:
            raise SetupError('Resolve the listed conflicts before setup.')
        print('stage=plan status=complete')
        if not args.apply:
            print('Dry run complete. App credentials and installation have not been verified.')
            return 0

        stage = 'settings'
        print('stage=settings status=started')
        rechecked = inspect_settings(gh, args.repo, args.reviewer)
        _unchanged(snapshot, rechecked)
        verified = apply_settings(gh, rechecked, plan, completed)
        _verify_settings(verified, app_id)
        print('stage=settings status=complete')

        stage = 'credentials'
        print('stage=credentials status=started')
        if app_service is None:
            import release_app_setup
            app_service = release_app_setup
        if args.create_app:
            # Registration can take minutes. Recheck before opening its browser.
            _unchanged(verified, inspect_settings(gh, args.repo, args.reviewer))
            registration = app_service.register_app(args.repo, verified['repository']['owner_type'])
            app_id = registration['id']
            pem = registration.pop('pem')
            completed.append('app:registered')
            print(f'App ID: {app_id}')
            print(f"Install: https://github.com/apps/{registration['slug']}/installations/new")
            print(f'Select {args.repo} only. Obtain organization approval if required.')
            print('Complete installation in GitHub, then type continue here.')
            stage = 'installation'
            print('stage=installation status=waiting')
            try:
                continuation = input('Type continue after installation: ')
            except EOFError:
                raise SetupError('Installation continuation was not received. Setup is incomplete.') from None
            if continuation.strip() != 'continue':
                raise SetupError('Installation continuation was not received. Setup is incomplete.')
        else:
            pem = _read_private_key(args.private_key_file)

        stage = 'app-verification'
        print('stage=app-verification status=started')
        if args.create_app:
            public_app = app_service.verify_installation(args.repo, app_id, pem, newly_created=True)
        else:
            public_app = app_service.verify_installation(args.repo, app_id, pem)
        for notice in public_app['notices']:
            print(f'Notice: {notice}')
        print('stage=app-verification status=complete')

        stage = 'credential-recheck'
        print('stage=credential-recheck status=started')
        current = inspect_settings(gh, args.repo, args.reviewer)
        _unchanged(verified, current)
        _verify_settings(current, app_id)
        print('stage=credential-recheck status=complete')
        if not current['secret_present']:
            stage = 'secret-upload'
            print('stage=secret-upload status=started')
            print('Use one setup operator at a time. GitHub has no atomic create-only secret upload.')
            upload_attempted = True
            gh.upload_key(args.repo, pem)
            uploaded = True
            completed.append('secret:RELEASE_APP_PRIVATE_KEY')
            print('stage=secret-upload status=complete')
        else:
            print('The stored secret was not verified. Verification covers only the supplied key.')
        pem = None

        stage = 'app-id'
        print('stage=app-id status=started')
        expected = dict(current, secret_present=True)
        current = inspect_settings(gh, args.repo, args.reviewer)
        _unchanged(expected, current)
        if 'RELEASE_APP_ID' not in current['variables']:
            gh.write('POST', f'/repos/{args.repo}/actions/variables',
                     {'name': 'RELEASE_APP_ID', 'value': str(app_id)})
            completed.append('variable:RELEASE_APP_ID')
        print('stage=app-id status=complete')

        stage = 'metadata'
        print('stage=metadata status=started')
        final = inspect_settings(gh, args.repo, args.reviewer)
        expected['variables'] = dict(expected['variables'], RELEASE_APP_ID=str(app_id))
        _unchanged(expected, final)
        _verify_settings(final, app_id)
        print('stage=metadata status=complete')
        for notice in final['readiness_notices']:
            print(f'Notice: {notice}')
        print('stage=report status=complete')
        print('Settings and supplied App access are verified. A successful protected workflow trial is still required.')
        print('New release enable variables are false. Existing values were preserved; activation is a separate administrator action.')
        return 0
    except (SetupError, OSError) as error:
        print(f'stage={stage} status=failed')
        print(str(error) if isinstance(error, SetupError)
              else 'A local or transport operation failed; no error details were retained.')
        print('Completed operations: ' + (', '.join(completed) if completed else 'none'))
        if app_id is not None and args.apply:
            print(f'App ID: {app_id}')
        if uploaded:
            print('Key upload completed. If RELEASE_APP_ID is absent, manual metadata repair is required: '
                  'confirm the App identity, then create that variable with the public App ID above. '
                  'Do not replace the stored secret.')
        elif upload_attempted:
            print('Upload failed or its result is uncertain. Inspect secret names before recovery. '
                  'Generate a replacement private key in GitHub App settings only if the secret is absent. '
                  'Save it privately and rerun with --app-id and --private-key-file. '
                  'If the secret exists without RELEASE_APP_ID, confirm the App identity and perform manual metadata repair.')
        elif 'app:registered' in completed:
            print('Keep the registered App. Complete its installation. Generate a replacement private key in GitHub App settings, '
                  'and rerun with --app-id and --private-key-file.')
        return 1
    finally:
        pem = None


if __name__ == '__main__':
    raise SystemExit(main())
