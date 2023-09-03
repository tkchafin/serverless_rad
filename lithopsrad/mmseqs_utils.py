import os
import sys
import re
import tempfile
import subprocess as sp

import lithopsrad.sequence as seq


def parse_mmseqs(tmp_hits, tmp_centroids):
	hits = dict()
	depths = dict()
	centroids = list()

	for record in read_mmseq_hits(tmp_hits):
		hit_spl = record[0].split(";size=")
		q_spl = record[1].split(";size=")

		# update hits record for centroid
		if hit_spl[0] not in hits:
			hits[hit_spl[0]] = list()
		hits[hit_spl[0]].append(q_spl[0])

		# update depth for centroid
		if hit_spl[0] not in depths:
			depths[hit_spl[0]] = int(hit_spl[1])
		depths[hit_spl[0]] += int(q_spl[1])

	# adjust fasta headers with new depths
	for header, sequence in seq.read_fasta(tmp_centroids):
		h_spl = header.split(";size=")
		if h_spl[0] in depths:
			header=str(h_spl[0]) + ";size=" + str(depths[h_spl[0]])
		centroids.append([header,sequence])

	return(hits, centroids)


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
			for key in hits:
				oline=str(key) + "\t" + ",".join(hits[key]) + "\n"
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
				if line[0] != line[1]:
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

