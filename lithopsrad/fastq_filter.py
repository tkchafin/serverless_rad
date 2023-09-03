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
        self.input_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_chunks"]))
        self.output_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_edits"]))

        # runtime params for edit step 
        self.minlen = self.runtime_config["edit"]["minlen"]
        self.truncqual = self.runtime_config["edit"]["truncqual"]
        self.maxns = self.runtime_config["edit"]["maxns"]
        self.maxee = self.runtime_config["edit"]["maxee"]
        self.maxee_rate = self.runtime_config["edit"]["maxee_rate"]

        # define map function 
        self._func = FASTQFilter._filter_fastq


    def _get_iterdata(self, obj):
        data = super()._get_iterdata(obj)
        data.update({
            "minlen": self.minlen,
            "truncqual": self.truncqual,
            "maxns": self.maxns,
            "maxee": self.maxee,
            "maxee_rate": self.maxee_rate
        })
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
        fout = base + ".edit"

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
        chunk_id = os.path.basename(os.path.splitext(fastq)[0])
        os.remove(fout)
        return {
            "chunk": chunk_id,
            #"filter_path": filter_path,
            "filtered_size": filtered_records
        }

