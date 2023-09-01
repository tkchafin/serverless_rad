import os 
import sys 
import tempfile
import subprocess as sp 
from lithops import FunctionExecutor

from lithopsrad.module import Module
import lithopsrad.sequence as seq
import lithopsrad.utils as utils

from lithops.storage import Storage

class FASTQFilter(Module):
    def __init__(self, lithops_config, runtime_config):
        super().__init__(lithops_config, runtime_config)
        self.setup()


    def setup(self):
        # Remote paths 
        self.run_path = utils.fix_dir_name(self.runtime_config["remote_paths"]["run_path"])
        self.fastq_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_path"]))
        self.fastq_chunks_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_chunks"]))
        self.fastq_edits_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_chunks"]))

        # runtime params for edit step 
        self.minlen = self.runtime_config["edit"]["minlen"]
        self.truncqual = self.runtime_config["edit"]["truncqual"]
        self.maxns = self.runtime_config["edit"]["maxns"]
        self.maxee = self.runtime_config["edit"]["maxee"]
        self.maxee_rate = self.runtime_config["edit"]["maxee_rate"]
    

    def run(self):

        chunks = self.list_remote_files(self.fastq_chunks_path)
        iterdata = [self._get_iterdata(chunk) for chunk in chunks]

        with FunctionExecutor(config=self.lithops_config) as fexec:
            # Passing the attributes as arguments
            fexec.map(FASTQFilter._filter_fastq, iterdata)
            results = fexec.get_result()
            print(results)


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
        """
        Filters fasta files using the vsearch tool.
        
        Args:
        - obj (CloudObject): The cloud object representing the chunk.
        - config (dict): The lithops configuration.
        - bucket (str): The bucket where the object is stored.
        - remote_path (str): The key/path of the object in the bucket.
        - minlen (int): Minimum length parameter for vsearch.
        - truncqual (int): Truncation quality parameter for vsearch.
        - maxns (int): Maximum Ns parameter for vsearch.
        - maxee (float): Maximum expected errors parameter for vsearch.
        - maxee_rate (float): Maximum expected error rate parameter for vsearch.
        - tmpdir (str, optional): Directory to use for temporary files. If not specified, uses the system default.
        """
        
        tmpdir = tmpdir or os.path.realpath(tempfile.gettempdir())
        os.chdir(tmpdir)

        fastq = obj.key  # Assuming CloudObject has a 'key' attribute for the object key/path
        tmp_path = os.path.join(tmpdir, os.path.basename(fastq))
        
        # Download the file to the temp directory
        utils._download_file(config, bucket, fastq, tmp_path)

        prefix = os.path.splitext(tmp_path)[0]
        base = os.path.basename(prefix)
        fout = prefix + ".temp.edit"

        cmd = [
            "vsearch",
            "-fastx_filter", tmp_path,
            "-fastqout", fout,
            "-fastq_minlen", str(minlen),
            "-fastq_truncqual", str(truncqual),
            "-fastq_maxns", str(maxns),
            "-fastq_maxee", str(maxee),
            "-fastq_maxee_rate", str(maxee_rate),
            "-threads", "1"  # You can modify this value as per your needs.
        ]

        proc = sp.Popen(cmd, stderr=sp.STDOUT, stdout=sp.PIPE, close_fds=True)
        res = proc.communicate()[0].decode("utf-8")
        print(res)

