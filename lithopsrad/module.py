import os 
import sys 
import time 
import pandas as pd 
from lithops import FunctionExecutor

import lithopsrad.utils as utils

class Module:
    def __init__(self, lithops_config, runtime_config):
        self.lithops_config = lithops_config
        self.runtime_config = runtime_config
        self.runtime = -1
        self._results = None
        self.bucket = self.runtime_config["global"]["bucket"]
        self.tmpdir = self.runtime_config["remote_paths"]["tmpdir"]
        self.nthreads = self.runtime_config["global"]["nthreads"]


    def validate(self):
        """
        Validate inputs
        """
        raise NotImplementedError("Subclasses should implement this method!")


    def run(self):
        """
        Method that implements the main functionality of the module.
        This should be overridden by subclasses.
        """
        raise NotImplementedError("Subclasses should implement this method!")


    @property
    def result(self):
        return self._get_result_as_df()


    def _get_result_as_df(self):
        # If results is None, return None
        if self._results is None:
            return None

        # If results is already a DataFrame, return it
        if isinstance(self._results, pd.DataFrame):
            return self._results

        # If results is a list, convert to DataFrame
        if isinstance(self._results, list):
            return pd.DataFrame(self._results)

        # If results is a dict, convert to DataFrame
        if isinstance(self._results, dict):
            return pd.DataFrame([self._results])

        # You can expand for other data types or add an error for unsupported types
        raise ValueError(f"Unsupported data type: {type(self._results)}")


    def list_remote_files(self, prefix=None):
        """
        Lists all the remote files (keys) in a specified bucket with an optional prefix.
        
        Args:
        - prefix (str, optional): Key prefix for filtering the listed files. If not provided,
                                all files in the bucket are listed.

        Returns:
        - list[str]: A list of object keys (file names) present in the specified bucket, filtered
                    by the provided prefix if any.
        """
        return utils._list_remote_files(self.lithops_config, self.bucket, prefix)


    def check_remote_files(self, prefix, subset=3):
        """
        Check if a subset of remote files under a given prefix is reachable.

        Args:
        - prefix (str): The prefix under which to check the files.
        - subset (int, optional): Number of files to check. Defaults to 3.

        Returns:
        - List of files that were checked.
        """
        return utils._check_remote_files(self.lithops_config, self.bucket, prefix, subset)


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
    

    @staticmethod
    def time_it(func):
        def wrapper(*args, **kwargs):
            start = time.time()
            result = func(*args, **kwargs)
            end = time.time()
            # Store the runtime in the instance (args[0] is 'self' for class methods)
            args[0].runtime = end - start
            return result
        return wrapper