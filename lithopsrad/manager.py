import os 
import sys

import json
import numpy as np
import pandas as pd 
import traceback

# import modules 
from lithopsrad.fastq_chunker import FASTQChunker
from lithopsrad.fastq_filter import FASTQFilter
from lithopsrad.fastq_derep import FASTQDerep
from lithopsrad.cluster_map import ClusterMap
from lithopsrad.cluster_merge import ClusterMerge

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
                print(traceback.format_exc())
                sys.exit()
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
        # fastq processing 
        self.run_fastq_chunker()
        self.run_fastq_filter()
        self.run_fastq_derep()

        # within-sample clustering 
        self.run_clust_within()
        self.run_clustmerge_within()

        # among-sample cluster merge
        self.run_clustmerge_across()
        
        # alignment
        # calling 
        # locus/catalog filter 

        # TODO: Implement cleanup() methods in relevant steps to remove intermediate files from bucket 


        sample_summary, chunk_summary = self.summarize_results()
        print(chunk_summary)
        print(sample_summary)


    @step_handler("FASTQChunker")
    def run_fastq_chunker(self):
        module = FASTQChunker(self.lithops_config, self.runtime_config)
        module.validate()
        return module
        

    @step_handler("FASTQFilter")
    def run_fastq_filter(self):
        module = FASTQFilter(self.lithops_config, self.runtime_config)
        module.validate()
        return module


    @step_handler("FASTQDerep")
    def run_fastq_derep(self):
        module = FASTQDerep(self.lithops_config, self.runtime_config)
        module.validate()
        return module
    
    @step_handler("ClusterMapWithin")
    def run_clust_within(self):
        module = ClusterMap(self.lithops_config, self.runtime_config, mode="clust_within")
        module.validate()
        return module

    @step_handler("ClusterMergeWithin")
    def run_clustmerge_within(self):
        module = ClusterMerge(self.lithops_config, self.runtime_config, mode="clust_within")
        module.validate()
        return module
    
    @step_handler("ClusterMergeAcross")
    def run_clustmerge_across(self):
        module = ClusterMerge(self.lithops_config, self.runtime_config, mode="clust_across")
        module.validate()
        return module


    def summarize_results(self):
        """
        Merges results into sample_summary and chunk_summary dataframes.
        Assumes the existence of self.results populated with dataframes.
        Returns: sample_summary and chunk_summary DataFrames
        """
        
        # Initialize the sample_summary and chunk_summary dataframes using the data from FASTQChunker.
        fastq_chunker_df = self.results["FASTQChunker"].copy()

        # chunk_summary initialization
        chunk_summary = fastq_chunker_df

        # sample_summary initialization
        sample_grouped = fastq_chunker_df.groupby('sample').agg({
            'chunk': 'count',
            'size': 'sum'
        }).rename(columns={'chunk': 'chunks'}).reset_index()
        sample_summary = sample_grouped[['sample', 'chunks', 'size']]

        # PLaceholder for catalog 
        if 'catalog' not in sample_summary['sample'].values:
            default_catalog = {'sample': 'catalog', 'chunks': np.nan, 'size': np.nan}
            sample_summary = pd.concat([sample_summary, pd.DataFrame([default_catalog])], ignore_index=True)

        # Iterate over each result in self.results
        for step, result in self.results.items():
            if step != "FASTQChunker":
                if 'sample' in result.columns:
                    sample_summary = pd.merge(sample_summary, result, on="sample", how="outer")
                if 'chunk' in result.columns:
                    chunk_summary = pd.merge(chunk_summary, result, on="chunk", how="outer")

        return sample_summary, chunk_summary


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

