import os 
import sys 

from lithops import FunctionExecutor
from lithops.storage import Storage

import lithopsrad.utils as utils

class Module:
    def __init__(self, lithops_config, runtime_config):
        self.lithops_config = lithops_config
        self.runtime_config = runtime_config
        self.bucket = self.runtime_config["global"]["bucket"]
    

    def run(self):
        """
        Method that implements the main functionality of the module.
        This should be overridden by subclasses.
        """
        raise NotImplementedError("Subclasses should implement this method!")
    

    def upload_file(self, remote_path, local_path, overwrite=True):
        """
        Upload a local file to a remote storage location.

        Args:
        - remote_path (str): The path in the remote storage to save the file.
        - local_path (str): The local path to the file to be uploaded.
        - overwrite (bool, optional): If True, overwrite the file if it exists in remote storage. Defaults to True.
        """
        
        print(f"Uploading file {local_path}... ", flush=True)

        # Check if the file exists in remote storage
        exists = utils._remote_file_exists(self.lithops_config, self.bucket, remote_path)

        # Handle file overwriting scenarios
        if exists and overwrite:
            print("File exists. Overwriting... ")
        elif exists and not overwrite:
            print("File exists. Skipping.")
            return utils._get_cloudobject(self.lithops_config, self.bucket, remote_path)

        # Upload the file
        return utils._upload_file(self.lithops_config, self.bucket, remote_path, local_path)


    def download_file(self, remote_path, local_path, overwrite=True):
        """
        Download a file from remote storage to a local location.

        Args:
        - remote_path (str): The path in the remote storage to fetch the file from.
        - local_path (str): The local path to save the downloaded file.
        - overwrite (bool, optional): If True, overwrite the file if it exists locally. Defaults to True.
        """
        
        print(f"Downloading file {remote_path}... ", end='', flush=True)

        # Check if the file exists locally
        if os.path.isfile(local_path):
            if not overwrite:
                print("File exists. Skipping.")
                return
            else:
                print("File exists. Overwriting... ", end='')

        # Fetch the file from remote storage and save locally
        utils._download_file(self.lithops_config, self.bucket, remote_path, local_path)

        print("Done.")


    def stream_file(self, remote_path):
        """
        Stream a file from remote storage line by line.

        Args:
        - remote_path (str): The path in the remote storage to stream the file from.
        
        Yields:
        - str: Next line from the file.
        """
        yield from utils._stream_file(self.lithops_config, self.bucket, remote_path)