"""
TODO: Multiple issues here (1) with large numbers of chunks
the file names can get really large; (2) mmseqs on localhost
can have conflicts trying to create intermediate databases, so
add a random identifier for each left-right pair here, which
will be used to name the outputs in the reducer 
"""
def binary_reducer_iterdata(queue):
	iterdata = list()
	kept = list()
	samples=dict()
	for item in queue:
		#print(item)
		if item["sample"] not in samples:
			samples[item["sample"]] = list()
		samples[item["sample"]].append(item)
	for sample in samples:
		while len(samples[sample]) > 1:
			iterdata.append({"left" : samples[sample].pop(0), "right" : samples[sample].pop(0)})
		if len(samples[sample]) > 0:
			kept.append(samples[sample].pop())
	return(kept, iterdata)

def make_merged_hits_table(left_hits, right_hits, infile, outfile):
	# get hit dicts for left, right, and joined outputs
	left_members = read_hits(left_hits)
	right_members = read_hits(right_hits)
	intermediate_hits = parse_hits(infile)

	joined_hits=dict()

	# grab constitutents for each centroid involved in a hit
	for int_hit in intermediate_hits:
		centroids=intermediate_hits[int_hit]
		joined_hits[int_hit] = centroids
		if int_hit in left_members:
			joined_hits[int_hit].extend(left_members[int_hit])
		if int_hit in right_members:
			joined_hits[int_hit].extend(right_members[int_hit])
		for c in centroids:
			if c in left_members:
				joined_hits[int_hit].extend(left_members[c])
			if c in right_members:
				joined_hits[int_hit].extend(right_members[c])
	#print(joined_hits)
	write_hits(joined_hits, outfile)
	return(joined_hits)