import boto3
import sys
from collections import namedtuple
import configparser
import os
import pathlib

from botocore.exceptions import ClientError

def get_credentials(profile):
    __s3_creds = namedtuple(
        "s3_credentials",
        ["access_key", "secret_key", "endpoint", "region", "profile_name"],)

    credential_file = configparser.ConfigParser()

    credential_file.read_file(open(os.path.expanduser("~/.aws/credentials"), "rt"))

    profile = "climb" if not profile else profile

    endpoint = "https://s3.climb.ac.uk"

    region = "s3"

    if credential_file:
        access_key = credential_file[profile]["aws_access_key_id"]
        secret_key = credential_file[profile]["aws_secret_access_key"]

    if os.getenv("AWS_ACCESS_KEY_ID"):
        access_key = os.getenv("AWS_ACCESS_KEY_ID")

    if os.getenv("AWS_SECRET_ACCESS_KEY"):
        secret_key = os.getenv("AWS_SECRET_ACCESS_KEY")

    if not access_key or not secret_key:
        error = """CLIMB S3 credentials could not be found, please provide valid credentials in one of the following ways:
            - In a correctly formatted config file (~/.aws/credentials)
            - As environmental variables 'AWS_ACCESS_KEY_ID' and 'AWS_SECRET_ACCESS_KEY'
            - As a command line argument, see --help for more details
        """
        print(error, file=sys.stderr)
        sys.exit(1)

    s3_credentials = __s3_creds(
        access_key=access_key,
        secret_key=secret_key,
        endpoint=endpoint,
        region=region,
        profile_name=profile,
    )

    return s3_credentials

def create_client(creds):
    s3_client = boto3.client("s3", endpoint_url=creds.endpoint,
    aws_access_key_id=creds.access_key,
    aws_secret_access_key=creds.secret_key)

    return s3_client

def _is_file_in_s3(client, folder, file):
    try:
        client.get_object(Bucket=folder, Key=file)
        return True
    except ClientError:
        return False


def is_file_in_s3(bucket, file, profile="climb"):
    creds = get_credentials(profile)
    client = create_client(creds)

    return _is_file_in_s3(client, bucket, file)
