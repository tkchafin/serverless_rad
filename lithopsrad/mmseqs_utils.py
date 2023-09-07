import os
import sys
import re
import tempfile
import subprocess as sp
import io
import re

import lithopsrad.sequence as seq


def get_cluster_info(fasta_file):
    """Get number of clusters and average depth from fasta file.
    
    Args:
        fasta_file (str): Path to the fasta file.
    
    Returns:
        int, float: Number of clusters and average depth.
    """
    
    num_clusters = 0
    total_depth = 0

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                num_clusters += 1
                size_value = line.split(";size=")[1]
                total_depth += int(size_value)

    average_depth = total_depth / num_clusters if num_clusters else 0

    return num_clusters, average_depth


def get_size_from_key(key):
    # Extracting the size value using regex
    match = re.search(r'size=(\d+)', key)
    if match:
        return int(match.group(1))
    return 0


# def make_merged_hits_table(left_hits, right_hits, infile):
#     # get hit dicts for left, right, and joined outputs
#     left_members = read_hits(left_hits)
#     right_members = read_hits(right_hits)
#     intermediate_hits = read_hits(infile)

#     joined_hits=dict()

#     # grab constitutents for each centroid involved in a hit
#     for int_hit in intermediate_hits:
#         centroids=intermediate_hits[int_hit]
#         joined_hits[int_hit] = centroids
#         if int_hit in left_members:
#             joined_hits[int_hit].extend(left_members[int_hit])
#         if int_hit in right_members:
#             joined_hits[int_hit].extend(right_members[int_hit])
#         for c in centroids:
#             if c in left_members:
#                 joined_hits[int_hit].extend(left_members[c])
#             if c in right_members:
#                 joined_hits[int_hit].extend(right_members[c])
#     return(joined_hits)


def make_merged_hits_table(left_hits, right_hits, infile):
    # get hit dicts for left, right, and joined outputs
    left_members = read_hits(left_hits)
    right_members = read_hits(right_hits)
    intermediate_hits = read_hits(infile)

    joined_hits = {}
    seen_centroids = set()

    # grab constituents for each centroid involved in a hit
    for int_hit, centroids in intermediate_hits.items():
        seen_centroids.add(int_hit)
        seen_centroids.update(centroids)

        for c in centroids:
            joined_hits.setdefault(int_hit, set()).update(left_members.get(c, []))
            joined_hits.setdefault(int_hit, set()).update(right_members.get(c, []))

        joined_hits.setdefault(int_hit, set()).update(left_members.get(int_hit, []))
        joined_hits.setdefault(int_hit, set()).update(right_members.get(int_hit, []))
        joined_hits[int_hit].update(centroids)

    # Add missing keys from left_members and right_members to joined_hits
    for key, value in left_members.items():
        if key not in seen_centroids:
            joined_hits.setdefault(key, set()).update(value)

    for key, value in right_members.items():
        if key not in seen_centroids:
            joined_hits.setdefault(key, set()).update(value)

    # Convert sets back to lists
    for key in joined_hits:
        joined_hits[key] = list(joined_hits[key])

    return joined_hits


def read_hits(infile):
    """
    Read hits table to a dict
    format: key \t member1,member2..memberN
    Args:
        infile(str): .temp.h file (hits file written by write_hits)
    Returns:
        Dict[str, List[str]]: Where key[str] is the derep ID for a centroid and List[str] are members
    Raises:
        FileNotFoundError: If hits file could not be found.
        IOError: If hits file could not be read from.
    """
    hits = {}
    
    try:
        with open(infile, "r") as ifh:
            for line in ifh:
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                split_line = line.split()

                # Handle only one element or empty second element in the line
                if len(split_line) == 1 or not split_line[1]:
                    hits[split_line[0]] = []
                else:
                    hits[split_line[0]] = split_line[1].split(",")
    except (IOError, FileNotFoundError):
        sys.exit(f"Could not open hits file {infile}")
    finally:
        return hits


def parse_mmseqs(tmp_hits, tmp_centroids):
    hits = dict()
    depths = dict()
    centroids = {}

    for record in read_mmseq_hits(tmp_hits):
        hit_spl = record[0].split(";size=")
        q_spl = record[1].split(";size=")

        if hit_spl[0] not in depths:
            depths[hit_spl[0]] = int(hit_spl[1])

        if hit_spl[0] == q_spl[0]:
            hits[hit_spl[0]] = list()
        else:
            # update hits record for centroid
            if hit_spl[0] not in hits:
                hits[hit_spl[0]] = list()
            hits[hit_spl[0]].append(q_spl[0])
            # update depth for centroid
            depths[hit_spl[0]] += int(q_spl[1])

    # adjust fasta headers with new depths
    for header, sequence in seq.read_fasta(tmp_centroids):
        h_spl = header.split(";size=")
        if h_spl[0] in depths:
            header = str(h_spl[0]) + ";size=" + str(depths[h_spl[0]])
        else:
            depths[h_spl[0]] = int(h_spl[1])

        centroids[header] = sequence

        # Ensure every centroid is represented in the hits table
        if h_spl[0] not in hits:
            hits[h_spl[0]] = []

    return hits, centroids


def write_hits(hits, outfile):
    """Write hits dictionary to file
    format: key \t member1,member2..memberN)

    Args:
        hits(Dict[str, List[str]]): Where key[str] is the derep ID for a centroid and List[str] are members
        outfile(str): Name for output hits file (usually .temp.h)
    Returns:
        None
    Raises:
        OSError: Unable to write to specified file
    """
    try:
        with open(outfile, "w") as ofh:
            # Sort hits dictionary by the length of their value lists in descending order
            sorted_hits = sorted(hits.items(), key=lambda x: len(x[1]), reverse=True)
            
            for key, value in sorted_hits:
                oline = str(key) + "\t" + ",".join(value) + "\n"
                ofh.write(oline)
    except OSError:
        sys.exit(f"Could not open file {outfile}")



def read_mmseq_hits(infile):
	try:
		with open(infile, "r") as ifh:
			for line in ifh:
				line = line.strip()
				# skip empty lines
				if not line:
					continue
				line=line.split()
				yield(line)
	except (IOError, FileNotFoundError):
		sys.exit(f"Could not open hits file {infile}")


def count_hits_from_stream(file_data):
    lines = file_data.strip().split("\n")
    total_entries = 0
    for line in lines:
        parts = line.split("\t")
        centroid = parts[0]
        hits = parts[1].split(",") if len(parts) > 1 else []
        total_entries += 1 + len(hits)  # 1 for centroid, rest for hits
    return total_entries / len(lines) if len(lines) else 0  # Average


def average_cluster_size(file_stream=None, file_path=None):
    if file_path is not None:
        with open(file_path, "r") as file:
            file_stream = file.read()
            return count_hits_from_stream(file_stream)
    else:
        return count_hits_from_stream(file_stream)

