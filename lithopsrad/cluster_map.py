import os 
import sys 
import tempfile
import subprocess as sp 
import shutil

from lithops import FunctionExecutor

from lithopsrad.module import Module
import lithopsrad.sequence as seq
import lithopsrad.utils as utils
import lithopsrad.mmseqs_utils as mmseqs_utils

class ClusterMap(Module):
    def __init__(self, lithops_config, runtime_config, mode="clust_within"):
        super().__init__(lithops_config, runtime_config)
        self.setup(mode)


    def setup(self, mode="clust_within"):
        # Remote paths 
        self.run_path = utils.fix_dir_name(self.runtime_config["remote_paths"]["run_path"])
        self.input_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["fastq_dereps"]))
        self.output_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["clust"]))

        # runtime params for cluster_map 
        self.cov = self.runtime_config[mode]["cov"]
        self.id = self.runtime_config[mode]["id"]
        self.cov_mode = self.runtime_config[mode]["cov_mode"]
        self.mask = self.runtime_config[mode]["mask"]
        self.mask_lower_case = self.runtime_config[mode]["mask-lower-case"]
        self.threads = self.runtime_config["global"]["nthreads"]

        # define function to run 
        self._func = ClusterMap._cluster_map


    def _get_iterdata(self, obj):
        data = super()._get_iterdata(obj)
        data.update({
            "cov": self.cov,
            "identity": self.id,
            "cov_mode": self.cov_mode,
            "mask": self.mask,
            "mask_lower_case": self.mask_lower_case,
            "threads": self.threads
        })
        return data

    @staticmethod
    def _cluster_map(obj, config, bucket, remote_path, cov, identity, cov_mode, mask, mask_lower_case, threads, tmpdir=None):
        # Set working directory
        tmpdir = tmpdir or os.path.realpath(tempfile.gettempdir())
        os.chdir(tmpdir)
        
        # Get the filename
        infile = obj['Key'] if isinstance(obj, dict) else obj.key
        tmp_path = os.path.join(tmpdir, os.path.basename(infile))
        
        # Download the file to the temp directory & configure paths 
        utils._download_file(config, bucket, infile, tmp_path)
        
        out_prefix = os.path.basename(os.path.splitext(tmp_path)[0])

        # Create a new sub-directory for MMSEQS2 temporary files
        mmseqs_tmp_dir = os.path.join(tmpdir, out_prefix)
        if os.path.exists(mmseqs_tmp_dir):
            shutil.rmtree(mmseqs_tmp_dir)  # Remove the directory if it exists
        os.makedirs(mmseqs_tmp_dir)
        
        # Run mmseqs
        cmd = [
            "mmseqs",
            "easy-linclust", 
            tmp_path,
            out_prefix,
            mmseqs_tmp_dir,
            "--min-seq-id", str(identity),
            "-c", str(cov),
            "--cov-mode", str(cov_mode),
            "--mask", str(int(mask)),
            "--createdb-mode", "0",
            "--add-self-matches", "0",
            "--mask-lower-case", str(int(mask_lower_case)),
            "--threads", str(threads), 
            "--remove-tmp-files", "1"
        ]
        
        proc = sp.Popen(cmd, stderr=sp.STDOUT, stdout=sp.PIPE, close_fds=True)
        res = proc.communicate()[0].decode("utf-8")
        print(" ".join(cmd))
        print(res)
        
        # Parse outputs into common hits-table format
        hits, centroids = mmseqs_utils.parse_mmseqs(out_prefix + "_cluster.tsv", out_prefix + "_rep_seq.fasta")
        os.remove(out_prefix + "_cluster.tsv")
        os.remove(out_prefix + "_rep_seq.fasta")

        # Write hits and centroids
        hits_path = os.path.join(mmseqs_tmp_dir, os.path.basename(out_prefix) + ".hits")
        centroids_path = os.path.join(mmseqs_tmp_dir, os.path.basename(out_prefix) + ".centroids")
        hout = os.path.join(remote_path, os.path.basename(out_prefix) + ".temp.hits")
        cout = os.path.join(remote_path, os.path.basename(out_prefix) + ".temp.centroids")
        
        # write files and grab results to report back 
        mmseqs_utils.write_hits(hits, hits_path)
        seq.write_fasta(centroids, centroids_path)
        centroids_num, cluster_depth = mmseqs_utils.get_cluster_info(centroids_path)

        # Upload the results
        utils._upload_file(config, bucket, hout, hits_path)
        utils._upload_file(config, bucket, cout, centroids_path)
        
        # Cleanup and return
        os.remove(hits_path)
        os.remove(centroids_path)
        shutil.rmtree(mmseqs_tmp_dir) 
        
        return {
            "chunk": out_prefix,
            "mean_depth_pre": cluster_depth,
            "clusters": centroids_num
        }