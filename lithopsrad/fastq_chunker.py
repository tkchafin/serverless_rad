import os 
import sys 

import pandas as pd 
from lithops import FunctionExecutor

from lithopsrad.module import Module, time_it
import lithopsrad.sequence as seq
import lithopsrad.utils as utils

class FASTQChunker(Module):
    def __init__(self, lithops_config, runtime_config):
        super().__init__(lithops_config, runtime_config)
        print(self.lithops_config)
        self.setup()


    def setup(self):
        self.input_fastq_dir = utils.fix_dir_name(self.runtime_config["input"]["input_fastq"])
        self.fastq_chunk_size = self.runtime_config["input"]["fastq_chunk_size"]
        self.ignore_R2 = self.runtime_config["input"]["ignore_R2"]

        self.run_path = utils.fix_dir_name(self.runtime_config["remote_paths"]["run_path"])
        self.fastq_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_path"]))
        self.output_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_chunks"]))

        self.overwrite_fastq = self.runtime_config["input"]["overwrite_fastq"]
        self.overwrite_chunks = self.runtime_config["input"]["overwrite_chunks"]

         # define map function 
        self._func = FASTQChunker._chunk_fastq
        self._reduce_func = FASTQChunker._chunk_fastq_reducer


    # overload Module validate 
    def validate(self):
        # check dir exists/ is readable 
        utils.check_file(self.input_fastq_dir)

        # get list of files, make sure isn't empty 
        local_files = [os.path.join(self.input_fastq_dir, file) for file in os.listdir(self.input_fastq_dir) if "R2" not in file]
        if not local_files:
            raise ValueError(f"No valid files found in the directory {self.input_fastq_dir}. Ensure that the files don't contain 'R2' in their names.")

        for file in local_files:
            # check if file exists/ is readable 
            utils.check_file(file)

            # check if it is possible fastq
            seq.check_fastq(file)

            # check if file is gzip 
            if utils.is_gzip(file):
                raise ValueError(f"File {file} appears to be a gzip compressed file, which is not supported.")
            

    # overload Module run 
    @time_it
    def run(self):
        # Upload local fastq files to bucket
        local_files = [os.path.join(self.input_fastq_dir, file) for file in os.listdir(self.input_fastq_dir) if "R2" not in file]
        cloud_paths = []
        for local_file in local_files:
            remote_path = os.path.join(self.fastq_path, os.path.basename(local_file))
            cloud_obj = self.upload_file(remote_path, local_file, overwrite=self.overwrite_fastq)
            cloud_paths.append(self._get_iterdata(obj=cloud_obj))

        # Chunk files w map_reduce              <-- (Currently ignores R2 reads)
        with FunctionExecutor(config=self.lithops_config) as fexec:
            fexec.map_reduce(self._func, 
                             cloud_paths, 
                             self._reduce_func,
                             chunksize=1,
                             obj_reduce_by_key=True,
                             obj_chunk_size=self.fastq_chunk_size, 
                             obj_newline="\n@")
            results = fexec.get_result()

            # check that record counts match after chunking 
            self._validate_chunks(results, local_files)
            self._results = results


    def _get_iterdata(self, obj):
        data = super()._get_iterdata(obj)
        data.update({
            "obj": obj,
        })
        return data


    def _validate_chunks(self, results, local_files):
        for local_file in local_files:
            # Obtain the base name of the local file
            base_name = os.path.basename(os.path.splitext(local_file)[0])

            # Count the FASTQ records in the local file
            local_records_count = seq.count_fastq_records(file_path=local_file)

            # Find the corresponding sample in the results
            matching_sample = next((res for res in results if res["sample"] == base_name), None)

            if matching_sample:
                if matching_sample["total_records"] != local_records_count:
                    raise ValueError(f"{base_name} has a mismatch in record counts. Local: {local_records_count}, Chunks: {matching_sample['total_records']}.")
            else:
                raise ValueError(f"{base_name} not found in results.")


    @staticmethod
    def _chunk_fastq_reducer(results):
        # assumed map_reduce was called with obj_reduce_by_key=True
        sample_id = os.path.basename(os.path.splitext(results[0]['original_key'])[0])
        chunk_paths = [item["chunk_path"] for item in results]
        chunks = [item["chunk"] for item in results]
        total_records = sum(res["record_count"] for res in results)
        sizes = [res["record_count"] for res in results]
        return {
            "sample": sample_id,
            "chunks": chunks,
            "chunk_paths": chunk_paths,
            "chunk_sizes" : sizes,
            "total_records": total_records
        }


    @staticmethod
    def _chunk_fastq(obj, config, bucket, remote_path, tmpdir=None):

        # Reading and counting the records
        data = obj.data_stream.read().decode('utf-8')
        record_count = seq.count_fastq_records(data)

        # Save chunk to remote storage 
        base_name = os.path.basename(obj.key)
        new_remote_path = os.path.join(remote_path, f"{obj.part}_{base_name}")
        chunk_cobj = utils._upload_file_from_stream(config, bucket, new_remote_path, data)

        chunk_id = os.path.basename(os.path.splitext(new_remote_path)[0])
        return {
            "original_key": os.path.basename(obj.key),
            "chunk_path": new_remote_path,
            "chunk": chunk_id,
            "record_count": record_count
        }


    def _get_result_as_df(self):
        # Check if _results is empty
        if not self._results:
            return None
        # Check if _results is already a DataFrame
        if isinstance(self._results, pd.DataFrame):
            return self._results

        all_data = []

        # Iterate over each dictionary in the list
        for item in self._results:
            sample = item['sample']

            chunks = item['chunks'] 

            # Fetch the sizes, and default to None if sizes are not provided.
            sizes = item.get('chunk_sizes', [None] * len(chunks))

            # Iterate over each chunk and its corresponding size
            for chunk, size in zip(chunks, sizes):
                all_data.append({
                    'sample': sample,
                    'chunk': chunk,
                    'size': size
                })

        # Transform the list of dictionaries into a DataFrame
        self._results = pd.DataFrame(all_data)
        return self._results


