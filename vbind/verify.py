import numpy as np
from tqdm import tqdm
from utils import ReverseComplement


def MatchesInFile(fname, search_string):
	# Find the number of lines in a file which contain a given string.
	occurrences = 0
	# print("Searching for a matching for {} in the file {}".format(search_string, fname))
	with open(fname, "r") as f:
		for line in f:
			if (line.strip(" ").strip("\n") == search_string):
				occurrences += 1
	return occurrences


def GetSubsequence(sequence, start, seqlen, is_circular):
	# Extract a subsequence from a string, from a start position, of a given length.
	# If the end position is less than the start, we will take a circular substring.
	subseq = ["N" for __ in range(seqlen)]
	for i in range(start, start + seqlen):
		if ((is_circular == 0) and (i >= len(sequence))):
			break
		subseq[i - start] = sequence[i % len(sequence)]
	return "".join(subseq)


def VerifyForwardMatching(rnap):
	# Verify the matching output from the dataset.
	# If the forward matching array element, F[i][j] = x, then we need to check if indeed the gene-substring gene[i:(i + lengths[i])] occurs x times in the pool.
	with open("./../data/input/%s" % (rnap.genefname), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	print("Genome: {}".format(gene_seq))

	for l in range(len(rnap.lengths)):
		is_valid = 1
		nuc_len = rnap.lengths[l][0]
		for i in tqdm(range(1, rnap.nForw.shape[1]), desc="Forward matchings for length %d" % (nuc_len)):
			nuc_len_idx = np.where(np.in1d(rnap.nForw[:, 0], rnap.lengths[l]))[0]
			gene_subseq = GetSubsequence(gene_seq, i - 1, nuc_len, rnap.is_circular)
			# gene_seq[(i - 1)  : (i - 1 + nuc_len) % len(gene_seq)]
			# print("Sub sequence from {} to {}: {}".format((i - 1), (i - 1 + nuc_len) % len(gene_seq), gene_subseq))
			given_matches = rnap.nForw[nuc_len_idx, i]
			found_matches = MatchesInFile("./../data/input/%s" % (rnap.poolfname), gene_subseq)
			
			if given_matches == found_matches:
				is_match = 1
				# print("[_/] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			else:
				is_match = 0
				print("[X] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			
			is_valid *= is_match
		
		if (is_valid == 1):
			print("vbind gave correct forward matching output for pool sequences with %d nucleotides." % (nuc_len))
		else:
			print("vbind gave incorrect forward matching output for pool sequences with %d nucleotides." % (nuc_len))
	return is_valid


def VerifyReverseMatching(rnap):
	# Verify the reverse matching output
	# If the reverse matching array element, R[i][j] = x, then we need to check if indeed the reverse complement of the gene-substring gene[i:(i + lengths[i])] occurs x times in the pool.
	with open("./../data/input/%s" % (rnap.genefname), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	for l in range(len(rnap.lengths)):
		is_valid = 1
		nuc_len = rnap.lengths[l][0]
		for i in tqdm(range(1, rnap.nRev.shape[1]), desc="Reverse matchings for length %d" % (nuc_len)):
			nuc_len_idx = np.where(np.in1d(rnap.nRev[:, 0], rnap.lengths[l]))[0]
			gene_subseq = ReverseComplement(GetSubsequence(gene_seq, i - 1, nuc_len, rnap.is_circular))
			given_matches = rnap.nRev[nuc_len_idx, i]
			found_matches = MatchesInFile("./../data/input/%s" % (rnap.poolfname), gene_subseq)
			
			if given_matches == found_matches:
				is_match = 1
				# print("[_/] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			else:
				is_match = 0
				print("[X] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			
			is_valid *= is_match
		
		if (is_valid == 1):
			print("vbind gave correct reverse matching output for pool sequences with %d nucleotides." % (nuc_len))
		else:
			print("vbind gave incorrect reverse matching output for pool sequences with %d nucleotides." % (nuc_len))
	return is_valid


def VerifyMatching(rnap):
	# Verify the forward and reverse matching output.
	is_valid = VerifyForwardMatching(rnap) * VerifyReverseMatching(rnap)
	return is_valid