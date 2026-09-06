"""Reuse exact digest-pinned images already present in the job's Docker daemon."""
import re
import subprocess


def ensure_pinned_image(image):
    if not re.fullmatch(r'ghcr\.io/[^\s@]+@sha256:[0-9a-f]{64}', image):
        raise ValueError('Image must use an immutable GHCR digest: ' + image)
    found = subprocess.run(['docker', 'image', 'inspect', image],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           check=False)
    if found.returncode == 0:
        print('stage=image_ready status=reused image=' + image, flush=True)
        return
    print('stage=image_ready status=pulling image=' + image, flush=True)
    subprocess.run(['docker', 'pull', image], check=True)
