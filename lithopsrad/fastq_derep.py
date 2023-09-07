import os 
import sys 
import tempfile
import subprocess as sp 
from lithops import FunctionExecutor

from lithopsrad.module import Module
import lithopsrad.sequence as seq
import lithopsrad.utils as utils

class FASTQDerep(Module):
    def __init__(self, lithops_config, runtime_config):
        super().__init__(lithops_config, runtime_config)
        self.setup()


    def setup(self):
        # Remote paths 
        self.run_path = utils.fix_dir_name(self.runtime_config["remote_paths"]["run_path"])
        self.input_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_edits"]))
        self.output_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_dereps"]))

        # runtime params for derep step 
        self.maxuniquesize = self.runtime_config["derep"]["maxuniquesize"]
        self.minuniquesize = self.runtime_config["derep"]["minuniquesize"]
        self.strand = self.runtime_config["derep"]["strand"]
        self.qmask = self.runtime_config["derep"]["qmask"]

        # define function to run 
        self._func = FASTQDerep._derep_fastq


    def _get_iterdata(self, obj):
        data = super()._get_iterdata(obj)
        data.update({
            "maxuniquesize": self.maxuniquesize,
            "minuniquesize": self.minuniquesize,
            "strand": self.strand,
            "qmask": self.qmask
        })
        return data


    @staticmethod
    def _derep_fastq(obj, config, bucket, remote_path, maxuniquesize, minuniquesize, strand, qmask, tmpdir=None):
        tmpdir = tmpdir or os.path.realpath(tempfile.gettempdir())
        os.chdir(tmpdir)

        fastq = obj['Key'] if isinstance(obj, dict) else obj.key
        tmp_path = os.path.join(tmpdir, os.path.basename(fastq))
        
        # 1. Download the file to the temp directory & configure paths 
        utils._download_file(config, bucket, fastq, tmp_path)

        prefix = os.path.splitext(tmp_path)[0]
        base = os.path.basename(prefix)
        fout = base + ".derep"
        label = base.split(".")[0] + "_d"

        # 2. Dereplicate using vsearch
        cmd = [
            "vsearch",
            "-fastx_uniques", tmp_path,
            "-fastqout", fout,
            "-strand", strand,
            "-relabel", label,
            "-sizeout",
            "-maxuniquesize", str(maxuniquesize),
            "-minuniquesize", str(minuniquesize),
            "-threads", "1"  # TODO: Change this if vthreads is available.
        ]
        proc = sp.Popen(cmd, stderr=sp.STDOUT, stdout=sp.PIPE, close_fds=True)
        res = proc.communicate()[0].decode("utf-8")
        print(res)

        # 3. Handle mask option and upload
        derep_path = os.path.join(remote_path, os.path.basename(fout))
        if qmask:
            mout = prefix + ".mask"
            cmd = [
                "vsearch",
                "-fastx_mask", fout,
                "-fastqout", mout,
                "-threads", "1"  # TODO: Change this if vthreads is available.
            ]
            proc = sp.Popen(cmd, stderr=sp.STDOUT, stdout=sp.PIPE, close_fds=True)
            res = proc.communicate()[0].decode("utf-8")
            print(res)

            # Upload masked file as .derep
            utils._upload_file(config, bucket, derep_path, mout)
            derep_size = seq.count_fastq_records(file_path=mout)
            os.remove(mout)  # remove the mask file after upload
        else:
            # Upload the dereplicated file
            utils._upload_file(config, bucket, derep_path, fout)
            derep_size = seq.count_fastq_records(file_path=fout)

        # 4. Clean up
        os.remove(fout)
        os.remove(tmp_path)

        # 5. Return the results
        return {
            "chunk": base,  
            "derep_size": derep_size
        }
