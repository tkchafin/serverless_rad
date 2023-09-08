import os 
import sys 
import tempfile
import subprocess as sp 
import shutil
import hashlib
from collections import defaultdict
from itertools import chain
from lithops import FunctionExecutor

from lithopsrad.module import Module
import lithopsrad.sequence as seq
import lithopsrad.utils as utils
import lithopsrad.mmseqs_utils as mmseqs_utils

class ClusterMerge(Module):
    def __init__(self, lithops_config, runtime_config, mode="clust_within"):
        super().__init__(lithops_config, runtime_config)
        self.setup(mode)


    def setup(self, mode="clust_within"):
        # Remote paths 
        self.run_path = utils.fix_dir_name(self.runtime_config["remote_paths"]["run_path"])
        self.input_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["clust"]))
        self.output_path = os.path.join(self.run_path, utils.fix_dir_name(self.runtime_config["remote_paths"]["clust"]))

        # runtime params for cluster_merge_pair
        self.cov = self.runtime_config[mode]["cov"]
        self.id = self.runtime_config[mode]["id"]
        self.cov_mode = self.runtime_config[mode]["cov_mode"]
        self.mask = self.runtime_config[mode]["mask"]
        self.mask_lower_case = self.runtime_config[mode]["mask-lower-case"]
        self.threads = self.runtime_config["global"]["nthreads"]

        # runtime params for cluster processing/ filtering 
        self.min_depth = self.runtime_config[mode]["min_depth"]
        self.max_depth = self.runtime_config[mode]["max_depth"]

        # define function to run 
        self._func = ClusterMerge._cluster_merge_pair
        self._process_func = ClusterMerge._process_clusters

        self.mode=mode


    def _get_iterdata(self, obj):
        data = super()._get_iterdata(obj)
        data.update({
            "cov": self.cov,
            "identity": self.id,
            "cov_mode": self.cov_mode,
            "mask": self.mask,
            "mask_lower_case": self.mask_lower_case,
            "threads": self.threads,
            "mode" : self.mode,
            "sample_file": True
        })
        chunk = os.path.basename(os.path.splitext(os.path.basename(obj))[0]).replace(".temp", "")
        try:
            chunk_id, sample = chunk.split("_", 1)
        except:
            chunk_id = chunk
            sample = chunk
        data.update({
            "chunk": chunk, 
            "chunk_id": chunk_id,
            "sample": sample
        })
        return data
    

    def _results_to_iterdata(self, results):
        """
        Formats the provided results into iterdata format and returns it.

        Args:
        - results (list[dict]): The results, typically a list of dictionaries.

        Returns:
        - list[dict]: A list of formatted iterdata.
        """
        formatted_iterdata = []
        for item in results:
            obj = item['hits_temp_path']
            data = super()._get_iterdata(obj)
            data.update({
                "cov": self.cov,
                "identity": self.id,
                "cov_mode": self.cov_mode,
                "mask": self.mask,
                "mask_lower_case": self.mask_lower_case,
                "threads": self.threads,
                "chunk": item['chunk'],
                "chunk_id": item['chunk_id'],
                "sample": item['sample'],
                "mode" : self.mode,
                "sample_file": False,
                "mean_depth_merged": item['mean_depth_merged'],
                "clusters_merged": item['clusters_merged'], 
                "hits_temp_path" : item["hits_temp_path"],
                "centroid_temp_path" : item["centroid_temp_path"]
            })
            formatted_iterdata.append(data)
        return formatted_iterdata


    def _get_process_iterdata(self, data):
        iterdata = []
        for item in data:
            it = {
                "config": self.lithops_config,
                "bucket": self.bucket,
                "remote_path": self.output_path,
                "tmpdir": self.tmpdir,
                'min_depth': self.min_depth,
                'max_depth': self.max_depth,
                'sample' : item["sample"]
            }
            if "hits_temp_path" in item:
                it.update({
                    'hits_temp_path':  item['hits_temp_path'],
                    'centroid_temp_path':  item['centroid_temp_path']
                })
            else:
                hits = str(utils._get_path(item["obj"]))
                centroids = hits.replace('.hits', '.centroids')
                it.update({
                    'hits_temp_path': hits,
                    'centroid_temp_path':  centroids
                })
            iterdata.append(it)
        return iterdata


    def run(self):
        # Check if _func is set
        if not self._func:
            raise NotImplementedError("Function to run not set for this module.")

        # get chunks to process 
        chunks = self.list_remote_files(self.input_path)
        if self.mode == "clust_within":
            # create iterdata 
            iterdata = [self._get_iterdata(chunk) for chunk in chunks if "temp" in str(chunk) and "hits" in str(chunk)]
            num_samples = len(set(data['sample'] for data in iterdata))
        else:
            # create iterdata 
            iterdata = [self._get_iterdata(chunk) for chunk in chunks if "temp" not in str(chunk) and "hits" in str(chunk)]
            for it in iterdata:
                it["sample"] = "catalog"
            num_samples = 1

        # Limit the iterdata to pairs for binary reduction and keep track of remaining pairs to be processed
        queue, pairs = self._binary_reducer_iterdata(iterdata)

        # run the function until all chunks reduced
        with FunctionExecutor(config=self.lithops_config) as fexec:
            while pairs:
                # generate unique id for filenames 
                for pair in pairs:
                    pair['pair_id'] = self._generate_filename(str(pair['left_obj']['chunk']) + str(pair['right_obj']['chunk']))

                fexec.map(self._func, pairs)
                results = fexec.get_result()

                # Add returned results to queue
                new_iterdata = self._results_to_iterdata(results)
                queue.extend(new_iterdata)

                # update queue 
                queue, pairs = self._binary_reducer_iterdata(queue)

                # if no pairs, exit loop
        
        # Check the number of result items
        if len(queue) != num_samples:
            raise Exception(f"Expected number of result items to be {num_samples}, but got {len(results)}")

        # Map process_cluster step 
        with FunctionExecutor(config=self.lithops_config) as fexec:
            process_iterdata = self._get_process_iterdata(queue)
            fexec.map(self._process_func, process_iterdata)
            self._results = fexec.get_result()


    def _generate_filename(self, input, length=15):
        """Generate a unique filename based on SHA-1 hashing."""
        hashed_name = hashlib.sha1(str(input).encode()).hexdigest()[:length]
        return f"{hashed_name}"

    # TODO: Would make better use of lambdas if we calculate sizes of each chunk 
    # and group reduce steps to contain as many chunks as possible, rather than binary reduce
    def _binary_reducer_iterdata(self, queue):
        print("BINARY_REDUCER")
        iterdata = []
        kept = []
        samples = {}
        for item in queue:
            if item["sample"] not in samples:
                samples[item["sample"]] = []
            samples[item["sample"]].append(item)
        print("Queue:")
        for key in samples.keys():
            print(key,":", len(samples[key]))
        for sample in samples:
            while len(samples[sample]) > 1:
                iterdata.append({
                    "left_obj": samples[sample].pop(0),
                    "right_obj": samples[sample].pop(0)
                })
            if len(samples[sample]) > 0:
                kept.append(samples[sample].pop())
        print("Pairs in work queue:", len(iterdata))
        print("Chunks remaining:", len(kept))
        print()
        return kept, iterdata

    @staticmethod 
    def _cluster_merge_pair(left_obj, right_obj, pair_id):

        # Extract main parameters from left_obj
        config = left_obj["config"]
        bucket = left_obj["bucket"]
        cov = left_obj["cov"]
        mode = left_obj["mode"]
        identity = left_obj["identity"]  # Ensure this key is consistent
        cov_mode = left_obj["cov_mode"]
        mask = left_obj["mask"]
        mask_lower_case = left_obj["mask_lower_case"]
        threads = left_obj["threads"]
        sample = left_obj["sample"]  
        remote_path = left_obj["remote_path"]
        tmpdir = left_obj["tmpdir"] or None

        # Set working directory
        tmpdir = tmpdir or os.path.realpath(tempfile.gettempdir())
        os.chdir(tmpdir)

        # Download files and set up local environment
        left_hits = str(utils._get_path(left_obj["obj"]))
        left_centroids = left_hits.replace('.hits', '.centroids')
        right_hits = str(utils._get_path(right_obj["obj"]))
        right_centroids = right_hits.replace('.hits', '.centroids')
        utils._download_file(config, bucket, left_centroids, os.path.join(tmpdir, str(pair_id)+"left.centroids"))
        utils._download_file(config, bucket, right_centroids, os.path.join(tmpdir, str(pair_id)+"right.centroids"))

        # Create a new sub-directory for MMSEQS2 temporary files
        out_prefix = os.path.join(tmpdir, str(pair_id))
        mmseqs_tmp_dir = os.path.join(tmpdir, out_prefix)
        if os.path.exists(mmseqs_tmp_dir):
            shutil.rmtree(mmseqs_tmp_dir) 
        os.makedirs(mmseqs_tmp_dir)

        # write concatenated centroids
        # TODO: Should we sort centroids before clustering?
        fasta_files = [os.path.join(tmpdir, str(pair_id)+"left.centroids"), os.path.join(tmpdir, str(pair_id)+"right.centroids")]
        joined_centroids = out_prefix+".joined.fasta"
        utils.concat_files(fasta_files, joined_centroids)
        os.remove(os.path.join(tmpdir, str(pair_id)+"left.centroids"))
        os.remove(os.path.join(tmpdir, str(pair_id)+"right.centroids"))

        # run clustering on joined centroids 
        cmd = [
            "mmseqs",
            "easy-linclust",
            joined_centroids,
            out_prefix,
            mmseqs_tmp_dir,
            "--min-seq-id", str(identity),
            "--add-self-matches", "0",
            "-c", str(cov),
            "--cov-mode", str(cov_mode),
            "--createdb-mode", "0",
            "--mask", str(int(mask)),
            "--mask-lower-case", str(int(mask_lower_case)),
            "--threads", str(threads), 
            "--remove-tmp-files", "1"
        ]
        proc = sp.Popen(cmd, stderr=sp.STDOUT, stdout=sp.PIPE, close_fds=True)
        res = proc.communicate()[0].decode("utf-8")
        print(" ".join(cmd))
        print(res)
        
        # Parse outputs into common hits-table format
        os.remove(joined_centroids)
        if mode == "clust_within":
            hits, centroids = mmseqs_utils.parse_mmseqs(out_prefix + "_cluster.tsv", 
                                                        out_prefix + "_rep_seq.fasta")
        else:
            # if clustering across, depth is calculated as number of centroids
            hits, centroids = mmseqs_utils.parse_mmseqs(out_prefix + "_cluster.tsv", 
                                                        out_prefix + "_rep_seq.fasta",
                                                        count_members = True)
        os.remove(out_prefix + "_cluster.tsv")
        os.remove(out_prefix + "_rep_seq.fasta")
        os.remove(out_prefix + "_all_seqs.fasta")

        # Upload hits results 
        int_hits_path = os.path.join(out_prefix + ".int.h")
        mmseqs_utils.write_hits(hits, int_hits_path)
        if mode == "clust_within":
            utils._download_file(config, bucket, left_hits, os.path.join(tmpdir, str(pair_id)+"left.hits"))
            utils._download_file(config, bucket, right_hits, os.path.join(tmpdir, str(pair_id)+"right.hits"))
        else:
            # if clustering across, ignore within-sample hits by creating empty hits files 
            if left_obj["sample_file"]:
                utils.touch_file(str(pair_id)+"left.hits")
            else:
                utils._download_file(config, bucket, left_hits, os.path.join(tmpdir, str(pair_id)+"left.hits"))
            if right_obj["sample_file"]:
                utils.touch_file(str(pair_id)+"right.hits")
            else:
                utils._download_file(config, bucket, right_hits, os.path.join(tmpdir, str(pair_id)+"right.hits"))
        joined_hits = mmseqs_utils.make_merged_hits_table(os.path.join(tmpdir, str(pair_id)+"left.hits"),
                                                          os.path.join(tmpdir, str(pair_id)+"right.hits"),
                                                          int_hits_path)
        os.remove(int_hits_path)
        os.remove(os.path.join(tmpdir, str(pair_id)+"left.hits"))
        os.remove(os.path.join(tmpdir, str(pair_id)+"right.hits"))
        hits_remote_path = os.path.join(remote_path, str(pair_id)+ ".hits")
        hits_temp_path = os.path.join(tmpdir, str(pair_id)+"joined.hits")
        mmseqs_utils.write_hits(joined_hits, hits_temp_path)
        utils._upload_file(config, 
                            bucket, 
                            hits_remote_path, 
                            os.path.join(tmpdir, str(pair_id)+"joined.hits"))
        os.remove(hits_temp_path)
        
        # upload centroids 
        centroids_temp_path = os.path.join(tmpdir, str(pair_id)+"joined.centroids")
        centroids_remote_path = os.path.join(remote_path, str(pair_id)+".centroids")
        seq.write_fasta(centroids, centroids_temp_path)
        utils._upload_file(config, 
                            bucket, 
                            centroids_remote_path, 
                            centroids_temp_path)
        centroids_num, cluster_depth = mmseqs_utils.get_cluster_info(centroids_temp_path)
        os.remove(centroids_temp_path)
        shutil.rmtree(mmseqs_tmp_dir) 

        # delete left and right files from buckets 
        if mode == "clust_within":
            for f in [left_hits, left_centroids, right_hits, right_centroids]:
                utils._delete_file(config, bucket, f)
        else:
            if not left_obj["sample_file"]:
                for f in [left_hits, left_centroids]:
                    utils._delete_file(config, bucket, f)
            if not right_obj["sample_file"]:
                for f in [right_hits, right_centroids]:
                    utils._delete_file(config, bucket, f)
            
        return {
            "chunk": pair_id,
            "chunk_id": pair_id,
            "sample": sample, 
            "centroid_temp_path": centroids_remote_path,
            "hits_temp_path": hits_remote_path,
            "mean_depth_merged": cluster_depth,
            "clusters_merged": centroids_num
        }
    
    def _process_clusters(config, bucket, remote_path, tmpdir, min_depth, max_depth, sample, hits_temp_path, centroid_temp_path):
        # Set working directory
        tmpdir = tmpdir or os.path.realpath(tempfile.gettempdir())
        os.chdir(tmpdir)

        # Dictionary to store valid FASTA records
        valid_records = {}
        current_header = None
        
        for line in utils._stream_file(config, bucket, centroid_temp_path):
            # If line is a header
            if line.startswith('>'):
                line = line.replace(">","")
                size = mmseqs_utils.get_size_from_key(line)
                
                # If size is within range
                if min_depth <= size <= max_depth:
                    current_header = line
                    valid_records[current_header] = ""
                else:
                    current_header = None
            elif current_header:  # If the line is part of a valid record
                valid_records[current_header] += line + '\n'

        # Writing valid FASTA records to a local file
        new_temp_file = os.path.join(tmpdir, f"{sample}_temp_centroids.fasta")
        seq.write_fasta(valid_records, new_temp_file)

        # Parse cluster number and sizes from the new file
        centroids_num, cluster_depth = mmseqs_utils.get_cluster_info(new_temp_file)

        # Upload new centroids file 
        centroids_new_path = os.path.join(remote_path, f"{sample}.centroids")
        utils._delete_file(config, bucket, centroid_temp_path)
        utils._upload_file(config, bucket, centroids_new_path, new_temp_file)
        os.remove(new_temp_file)

        # Rename the old hits file
        hits_new_path = os.path.join(remote_path, f"{sample}.hits")
        utils._rename_file(config, bucket, hits_temp_path, hits_new_path)

        # Format the results
        formatted_results = {
            "sample": sample,
            "mean_depth_merged": cluster_depth,
            "clusters_merged": centroids_num
        }
        
        return formatted_results

