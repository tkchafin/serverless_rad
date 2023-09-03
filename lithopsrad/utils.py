
import os 
from pathlib import Path

from lithops.storage import Storage
from lithops.storage.utils import CloudObject


def _check_remote_files(config, bucket, prefix, subset=None):
    """
    Check if a subset of remote files under a given prefix is reachable using the Lithops Storage API.
    
    Args:
    - config (dict): The lithops configuration.
    - bucket (str): The bucket where the files reside.
    - prefix (str): The prefix under which to check the files.
    - subset (int, optional): Number of files to check. Defaults to 3.
    
    Returns:
    - List of files that were checked.
    
    Raises:
    - ValueError: If any of the subset files is not accessible.
    """
    storage = Storage(config=config)
    files = storage.list_objects(bucket, prefix)
    
    # Get a subset of files
    if subset is not None:
        subset_files = files[:subset]
    else:
        subset_files = files
    
    for file in subset_files:
        try:
            # Checking file metadata is a lightweight way to see if it's reachable.
            storage.head_object(bucket, file['Key'])
        except Exception as e:
            raise ValueError(f"Unable to access file {file['Key']}. Reason: {str(e)}")
    
    return subset_files



def check_bucket(config, bucket):
    """
    Check if a bucket exists and is readable using the Lithops Storage API.
    
    Args:
    - config (dict): The lithops configuration.
    - bucket (str): The bucket to check.
    
    Raises:
    - ValueError: If the bucket does not exist.
    - PermissionError: If the bucket is not readable.
    """
    storage = Storage(config=config)
    
    try:
        # This is a lightweight operation to check if the bucket is accessible.
        storage.head_bucket(bucket)
    except Exception as e:
        raise PermissionError(f"Unable to access bucket {bucket}. Reason: {str(e)}")


def is_gzip(file_path):
    with open(file_path, 'rb') as f:
        file_start = f.read(2)
        return file_start == b'\x1f\x8b'


def check_file(file_path):
    """
    Checks if the given file exists, is readable, and conforms to FASTQ format.

    Args:
    - file_path (str): Path to the file to be checked.

    Raises:
    - ValueError: If the file doesn't exist.
    - PermissionError: If the file is not readable.
    """
    # Check file existence
    if not os.path.exists(file_path):
        raise ValueError(f"File {file_path} does not exist.")
    
    # Check file permissions
    if not os.access(file_path, os.R_OK):
        raise PermissionError(f"File {file_path} is not readable.")


def _list_remote_files(config, bucket, prefix=None):
    storage = Storage(config=config)        
    return storage.list_keys(bucket, prefix=prefix)


def _cloudobject_url(cobj):
    path = f'{cobj.backend}://{cobj.bucket}/{cobj.key}'
    return path


def _remote_file_exists(config, bucket, remote_path):
    """Check if the file exists in the specified bucket using lithops storage."""
    storage = Storage(config=config)
    try:
        storage.head_object(bucket, remote_path)
        return True
    except:
        return False


def _get_cloudobject(config, bucket, remote_path):
    """
    Construct a CloudObject for the given remote path in the specified bucket.
    
    Args:
    - config (dict): The lithops configuration.
    - bucket (str): The bucket where the object is stored.
    - remote_path (str): The key/path of the object in the bucket.
    
    Returns:
    - CloudObject: The constructed cloud object.
    """
    backend = Storage(config=config).backend
    return CloudObject(backend, bucket, remote_path)


def _upload_file_from_stream(config, bucket, remote_path, stream):
    """Upload a file stream to the specified bucket using lithops storage."""
    storage = Storage(config=config)
    storage.put_object(bucket, f'{remote_path}', stream)
    return _get_cloudobject(config, bucket, remote_path)


def _upload_file(config, bucket, remote_path, local_path):
    """Upload a file to the specified bucket using lithops storage."""
    storage = Storage(config=config)
    key = os.path.basename(local_path)
    with open(f'{local_path}', 'rb') as fl:
        storage.put_object(bucket, f'{remote_path}', fl)
    return _get_cloudobject(config, bucket, remote_path)


def _download_file(config, bucket, remote_path, local_path):
    """Download a file from the specified bucket using lithops storage."""
    storage = Storage(config=config)
    fobj = storage.get_object(bucket, remote_path)
    with open(local_path, "wb") as f:  # note the change from 'w' to 'wb' as we're writing bytes
        f.write(fobj)


def _stream_file(config, bucket, remote_path):
    """
    Generator function to iterate over file from remote storage
    
    Args:
    - config (dict): Lithops configuration.
    - bucket (str): The storage bucket name.
    - remote_path (str): The path in the remote storage to stream the file from.
    
    Yields:
    - str: Next line from the file.
    """
    storage = Storage(config=config)
    fobj = storage.get_object(bucket, remote_path)
    
    # Stream the file line by line
    for line in fobj.decode('UTF-8').splitlines():
        line = line.strip()
        if line:  # Ensure the line is not empty
            yield line


def fix_dir_name(d):
    """
    Ensure the directory path ends with a separator.
    """
    return str(Path(d).as_posix()) + '/'

