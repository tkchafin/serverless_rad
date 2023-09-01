import os 
import sys
import json

import lithopsrad.utils as utils
from lithopsrad.fastq_chunker import FASTQChunker


def step_error_handler(step_name):
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                print(f"Pipeline execution failed during {step_name}: {str(e)}")
        return wrapper
    return decorator


class PipelineManager:
    def __init__(self, config_file="lithops_config"):
        self.config_file = config_file
        self.lithops_config, self.runtime_config = self.get_params_from_json()


    def run(self):
        """
        Execute the pipeline.
        """
        self.run_fastq_chunker()


    @step_error_handler("FASTQChunker")
    def run_fastq_chunker(self):
        fastq_chunker = FASTQChunker(self.lithops_config, self.runtime_config)
        fastq_chunker.run()


    def get_params_from_json(self):
        """
        Read the configuration from a JSON file.
        """
        with open(self.config_file) as jfh:
            params = json.load(jfh)
        
        config = params["lithops_config"]
        args = params["runtime_args"]

        args = self.set_defaults(args)
        return config, args


    def set_defaults(self, args):
        """
        Set default parameters if not present in the config.
        """
        # set defaults
        pass
        return args



