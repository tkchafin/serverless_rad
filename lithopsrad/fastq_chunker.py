import os 
import sys 

from lithops import FunctionExecutor

from lithopsrad.module import Module
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
        self.fastq_chunks_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_chunks"]))

        self.overwrite_fastq = self.runtime_config["input"]["overwrite_fastq"]
        self.overwrite_chunks = self.runtime_config["input"]["overwrite_chunks"]
    

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
            fexec.map_reduce(FASTQChunker._chunk_fastq, 
                             cloud_paths, 
                             FASTQChunker._chunk_fastq_reducer,
                             chunksize=1,
                             obj_reduce_by_key=True,
                             obj_chunk_size=self.fastq_chunk_size, 
                             obj_newline="\n@")
            results = fexec.get_result()

            # check that record counts match after chunking 
            self._validate_chunks(results, local_files)
    

    def _validate_chunks(self, results, local_files):
        for local_file in local_files:
            # Obtain the base name of the local file
            base_name = os.path.basename(local_file)

            # Count the FASTQ records in the local file
            local_records_count = seq.count_fastq_records(file_path=local_file)

            # Find the corresponding sample in the results
            matching_sample = next((res for res in results if res["sample"] == base_name), None)

            if matching_sample:
                if matching_sample["total_records"] != local_records_count:
                    print(f"Error: {base_name} has a mismatch in record counts. Local: {local_records_count}, Chunks: {matching_sample['total_records']}.")
            else:
                print(f"Error: {base_name} not found in results.")


    def _get_iterdata(self, obj):
        data = {
            "obj": obj,
            "config": self.lithops_config,
            "bucket": self.bucket,
            "remote_path": self.fastq_chunks_path
        }
        return data


    @staticmethod
    def _chunk_fastq_reducer(results):
        # assumed map_reduce was called with obj_reduce_by_key=True
        original_key = results[0]['original_key']
        chunk_paths = [item["chunk_path"] for item in results]
        chunks = [item["chunk"] for item in results]
        total_records = sum(res["record_count"] for res in results)
        return {
            "sample": original_key,
            "chunks": chunks,
            "chunk_paths": chunk_paths,
            "total_records": total_records
        }



    @staticmethod
    def _chunk_fastq(obj, config, bucket, remote_path):

        # Reading and counting the records
        data = obj.data_stream.read().decode('utf-8')
        record_count = seq.count_fastq_records(data)

        # Save chunk to remote storage 
        base_name = os.path.basename(obj.key)
        new_remote_path = os.path.join(remote_path, f"{obj.part}_{base_name}")
        chunk_cobj = utils._upload_file_from_stream(config, bucket, new_remote_path, data)

        return {
            "original_key": os.path.basename(obj.key),
            "chunk_path": new_remote_path,
            "chunk": os.path.basename(new_remote_path),
            "record_count": record_count
        }



