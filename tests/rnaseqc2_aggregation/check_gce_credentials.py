"""Check the image selects VM credentials; replace only metadata-server I/O."""

from unittest.mock import patch

from gslib import gcs_json_credentials

credentials_type = gcs_json_credentials.credentials_lib.GceAssertionCredentials
with (
    patch.object(
        gcs_json_credentials, "GetGceCredentialCacheFilename", return_value=None,
    ),
    patch.object(
        credentials_type,
        "_ScopesFromMetadataServer",
        return_value=["https://www.googleapis.com/auth/devstorage.read_only"],
    ),
):
    credentials = gcs_json_credentials._GetGceCreds()
    assert isinstance(credentials, credentials_type), (
        "The container must select VM service-account credentials without a key file"
    )
print("VM service-account credential selection passed (metadata I/O replaced).")
