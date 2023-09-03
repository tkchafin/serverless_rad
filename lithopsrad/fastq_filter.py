import os 
import sys 
import tempfile
import subprocess as sp 
from lithops import FunctionExecutor

from lithopsrad.module import Module
import lithopsrad.sequence as seq
import lithopsrad.utils as utils

class FASTQFilter(Module):
    def __init__(self, lithops_config, runtime_config):
        super().__init__(lithops_config, runtime_config)
        self.setup()


    def setup(self):
        # Remote paths 
        self.run_path = utils.fix_dir_name(self.runtime_config["remote_paths"]["run_path"])
        self.fastq_chunks_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_chunks"]))
        self.fastq_edits_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_chunks"]))

        # runtime params for edit step 
        self.minlen = self.runtime_config["edit"]["minlen"]
        self.truncqual = self.runtime_config["edit"]["truncqual"]
        self.maxns = self.runtime_config["edit"]["maxns"]
        self.maxee = self.runtime_config["edit"]["maxee"]
        self.maxee_rate = self.runtime_config["edit"]["maxee_rate"]
    

    def validate(self):
        # Check if the bucket in which the chunks reside exists and is accessible.
        utils.check_bucket(self.lithops_config, self.bucket)

        # Check if the list of remote files is non-empty.
        chunks = self.list_remote_files(self.fastq_chunks_path)
        if not chunks:
            raise ValueError(f"No files found in {self.fastq_chunks_path}")

        # check remote files are accessible (only checking a subset of them)
        self.check_remote_files(self.fastq_chunks_path, subset=3)


    @Module.time_it
    def run(self):
        # get chunks to process 
        chunks = self.list_remote_files(self.fastq_chunks_path)

        # create iterdata 
        iterdata = [self._get_iterdata(chunk) for chunk in chunks]

        # run filter on each chunk
        with FunctionExecutor(config=self.lithops_config) as fexec:
            # Passing the attributes as arguments
            fexec.map(FASTQFilter._filter_fastq, iterdata)
            results = fexec.get_result()
            # TODO: Parse filtered chunks and return back sizes for logging 


    def _get_iterdata(self, obj):
        cloud_path = utils._get_cloudobject(self.lithops_config, self.bucket, obj)
        data = {
            "obj": cloud_path,
            "config": self.lithops_config,
            "bucket": self.bucket,
            "remote_path": self.fastq_edits_path,
            "minlen": self.minlen,
            "truncqual": self.truncqual,
            "maxns": self.maxns,
            "maxee": self.maxee,
            "maxee_rate": self.maxee_rate,
            "tmpdir": self.tmpdir
        }
        return data


    @staticmethod
    def _filter_fastq(obj, config, bucket, remote_path, minlen, truncqual, maxns, maxee, maxee_rate, tmpdir=None):
        tmpdir = tmpdir or os.path.realpath(tempfile.gettempdir())
        os.chdir(tmpdir)

        fastq = obj['Key'] if isinstance(obj, dict) else obj.key
        tmp_path = os.path.join(tmpdir, os.path.basename(fastq))
        
        # 1. Download the file to the temp directory & configure paths 
        utils._download_file(config, bucket, fastq, tmp_path)

        prefix = os.path.splitext(tmp_path)[0]
        base = os.path.basename(prefix)
        fout = prefix + ".temp.edit"

        # 2. Filter using vsearch
        cmd = [
            "vsearch",
            "-fastx_filter", tmp_path,
            "-fastqout", fout,
            "-fastq_minlen", str(minlen),
            "-fastq_truncqual", str(truncqual),
            "-fastq_maxns", str(maxns),
            "-fastq_maxee", str(maxee),
            "-fastq_maxee_rate", str(maxee_rate),
            "-threads", "1"
        ]

        proc = sp.Popen(cmd, stderr=sp.STDOUT, stdout=sp.PIPE, close_fds=True)
        res = proc.communicate()[0].decode("utf-8")
        print(res)

        # 3. Delete temp files and upload result to bucket 
        os.remove(tmp_path)

        filter_path = os.path.join(remote_path, os.path.basename(fout))
        utils._upload_file(config, bucket, filter_path, fout)

        # 4. Get the count of filtered FASTQ records
        filtered_records = seq.count_fastq_records(file_path=fout)

        # 5. Clean up and return results 
        os.remove(fout)
        return {
            "chunk_path": fastq,
            "filter_path": filter_path,
            "filtered_records": filtered_records
        }

