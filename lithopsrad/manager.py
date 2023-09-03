import os 
import sys

import json
import pandas as pd 

# import modules 
from lithopsrad.fastq_chunker import FASTQChunker
from lithopsrad.fastq_filter import FASTQFilter
from lithopsrad.fastq_derep import FASTQDerep


def step_handler(step_name):
    def decorator(func):
        def wrapper(self, *args, **kwargs): 
            try:
                # Obtain the module instance from the decorated function.
                module = func(self, *args, **kwargs)

                # Run the module.
                module.run()

                # Store the result in the results dictionary.
                self.results[step_name] = module.result

                # Return the module instance.
                return module

            except Exception as e:
                print(f"Pipeline execution failed during {step_name}: {str(e)}")
        return wrapper
    return decorator


class PipelineManager:
    def __init__(self, config_file="lithops_config"):
        self.config_file = config_file
        self.lithops_config, self.runtime_config = self.get_params_from_json()
        self.results = {}  # This will store the results of each module.


    def run(self):
        """
        Execute the pipeline.
        """
        self.run_fastq_chunker()
        self.run_fastq_filter()
        self.run_fastq_derep()


        merged_df = self.results["FASTQChunker"]
        for step, result in self.results.items():
            if step != "FASTQChunker":  # as we've already initialized with this
                merged_df = pd.merge(merged_df, result, on="chunk", how="outer")

        print(merged_df)


    @step_handler("FASTQChunker")
    def run_fastq_chunker(self):
        # Instantiate and validate.
        module = FASTQChunker(self.lithops_config, self.runtime_config)
        module.validate()
        return module
        

    @step_handler("FASTQFilter")
    def run_fastq_filter(self):
        # Instantiate and validate.
        module = FASTQFilter(self.lithops_config, self.runtime_config)
        module.validate()
        return module


    @step_handler("FASTQDerep")
    def run_fastq_derep(self):
        # Instantiate and validate.
        module = FASTQDerep(self.lithops_config, self.runtime_config)
        module.validate()
        return module


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


    def validate_args(self, args):
        pass


    def set_defaults(self, args):
        """
        Set default parameters if not present in the config.
        """
        # set defaults
        if "tmpdir" not in args["remote_paths"]:
            args["remote_paths"]["tmpdir"] = None
        if "nthreads" not in args["global"]:
            args["global"]["nthreads"] = 1
        return args

