import numpy as np
from utils import WriteMatrixToFile, IntArrayToString

def RecordNormalizedData(rnap):
	# Save the normalized data to a text file.
	million = 1E6
	alignments = ["linear", "circular"]
	is_nan = 0
	with open("./../data/output/normalized_%s_%s_%s_%d.txt" % (rnap.genefname, rnap.poolfname, alignments[rnap.is_circular], rnap.tolerance), "w") as fp:
		fp.write("# Normalized data for \n")
		fp.write("# Genome: %s\n# Pool: %s\n# Topology: %s\n# Tolerance: %d\n#\n" % (rnap.genefname, rnap.poolfname, alignments[rnap.is_circular], rnap.tolerance))	
		fp.write("-----------------------------------------------------------------------\n")
		fp.write("{:^6} | {:^8} | {:^8} | {:^12} | {:^8} | {:^12}\n".format("Length", "In pool", "Forward", "N(Forward)", "Reverse", "N(Reverse)"))
		fp.write("=======================================================================\n")
		for l in range(len(rnap.lengths)):
			total_seqs = np.sum([rnap.nSequencesByLengths[l] for l in rnap.lengths[l]])
			
			nuclens = np.where(np.in1d(rnap.nForw[:, 0], rnap.lengths[l]))[0]
			nForw_matchings = np.sum(rnap.nForw[nuclens, 1:])
			if (total_seqs > 0):
				norm_forw_matchings = nForw_matchings * million/total_seqs
			else:
				is_nan = 1
				norm_forw_matchings = -1
			
			nuclens = np.where(np.in1d(rnap.nRev[:, 0], rnap.lengths[l]))[0]
			nRev_matchings = np.sum(rnap.nRev[nuclens, 1:])
			if (total_seqs > 0):
				norm_rev_matchings = nRev_matchings * million/total_seqs
			else:
				norm_rev_matchings = -1
			
			fp.write("{:^6} | {:^8} | {:^8} | {:^12} | {:^8} | {:^12}\n".format(IntArrayToString(rnap.lengths[l]), total_seqs, "%d" % nForw_matchings, "%.2g" % norm_forw_matchings, "%d" % nRev_matchings, "%.2g" % norm_rev_matchings))
			fp.write("-----------------------------------------------------------------------\n")

		fp.write("#\n# N(.) refers to the normalized read, i.e., reads per million nucleotides.\n")
		if (is_nan == 1):
			fp.write("# -1: Not defined since the total number of sequences in the pool of the corresponding length is 0.")
	return None

def Normalize(rnap, forward, reverse):
	# Multiply the number of reads at every position by a constant that depends on the nucleatide length.
	# This is to ensure that all reads are scaled for a pool size that is 1,000,000.
	# Scaling array provides the scaling constants by pool gene tol length.
	million = 1E6
	rnap.normalized_forward = np.zeros((len(rnap.lengths), rnap.gene_length), dtype = np.float)
	rnap.normalized_reverse = np.zeros((len(rnap.lengths), rnap.gene_length), dtype = np.float)
	for s in range(len(rnap.lengths)):
		total_seqs = np.sum([rnap.nSequencesByLengths[l] for l in rnap.lengths[s]])
		rnap.normalized_forward[s, :] = forward[s, :] * million/total_seqs
		rnap.normalized_reverse[s, :] = reverse[s, :] * million/total_seqs
	WriteMatrixToFile(rnap.norm_ForwDataFile, rnap.normalized_forward, is_append = 0, dataType='f')
	WriteMatrixToFile(rnap.norm_RevDataFile, rnap.normalized_reverse, is_append = 0, dataType='f')
	return None