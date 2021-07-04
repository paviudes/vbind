import sys
from multiprocessing import Process, Array
import os
import time
import numpy as np
from startup import DisplayLogoLicense, CheckDependencies, Usage

DEBUG = 0
QUIET = 0

class RNAProfiler():
	"""
	RNAProfiler class: Contains all the information required to perform a round of Pattern matching
	"""

	def __init__(self):
		# Initialize all the variables in the scope of the RNAProfiler class
		self.interfaceFile = "NAN"
		self.poolfname = "NAN"
		self.genefname = "NAN"
		self.matchFrequencyFname = "NAN"
		self.dateStamp = ""

		self.rawForwMatchFname = ""
		self.rawRevMatchFname = ""

		self.forwardMatchFname = ""
		self.reverseMatchFname = ""

		self.runtime = 0

		self.gene_length = 0
		self.MAXCONCURRENT = 100000
		self.poolsize = 0
		self.tolerance = 0
		self.max_nuc_length = 0
		self.is_circular = 1
		self.file_skip = 4

		self.gene_matrix = np.array([])
		self.forwardMatches = np.array([])
		self.reverseMatches = np.array([])
		self.forwardMatchCounts = np.array([])
		self.reverseMatchCounts = np.array([])
		return None


# Load a RNAProfiler object
def Load(bindObj, interfaceFile = "NAN"):
	# Load all the information required in the RNAProfiler object for Pattern matching
	# Read the interface file and load the parameters for Binding
	if os.path.isfile(interfaceFile):
		bindObj.interfaceFile = interfaceFile
		if (QUIET == 0):
			print("\033[2mLoading information required for pattern matching from %s.\033[0m" % bindObj.interfaceFile)

		with open(bindObj.interfaceFile, 'r') as uiFid:
			# 1. file name of pstvd sequence
			# 2. file name of nucleotides in the pool
			# 3. tolerance
			# 4. circular or linear mapping
			bindObj.genefname = uiFid.readline().strip("\n").strip(" ")
			bindObj.poolfname = uiFid.readline().strip("\n").strip(" ")
			bindObj.tolerance = int(uiFid.readline().strip("\n").strip(" "))
			bindObj.is_circular = int(uiFid.readline().strip("\n").strip(" "))
	else:
		LoadRNAProfilerFromConsole(bindObj)
	return None


def LoadRNAProfilerFromConsole(bindObj):
	# Load a RNAProfiler object from taking console inputs
	print("\033[93mEnter the name of the file containing the sRNA pool.\033[0m")
	print(">>"),
	bindObj.poolfname = raw_input().strip("\n")

	print("\033[93mEnter the name of the file containing the PSTVd sequence.\033[0m")
	print(">>"),
	bindObj.genefname = raw_input().strip("\n")

	print("\033[93mEnter the maximum number of mismatches allowed.\033[0m")
	print(">>"),
	bindObj.tolerance = int(raw_input().strip("\n"))

	return None


def NucleotideEncoder(nucleotide, blockSize, reverse):
	# Convert a nucleotide sequence to a binary sequence using an alphabet mapping of the following kind.
	# A --> 1000, T --> 0100, G --> 0010, C --> 0001
	encoding = {'A':'1000', 'T':'0100', 'G':'0010', 'C':'0001'}
	reverseEncoding = {'A':'0100', 'T':'1000', 'G': '0001', 'C':'0010'}
	nucLen = len(nucleotide)
	binaryNucleotide = np.zeros(4*blockSize, dtype = int)
	if reverse == 0:
		for i in range(nucLen):
			binaryNucleotide[4*i : 4*(i + 1)] = list(map(int, encoding[nucleotide[i]]))
	else:
		for i in range(nucLen):
			binaryNucleotide[4*((nucLen - i) - 1): 4 * (nucLen - i)] = list(map(int, reverseEncoding[nucleotide[i]]))
	if DEBUG > 1:
		print("\033[2mNucleotide sequence: %s and encoding produced for direction %d:" % (nucleotide, reverse))
		print(binaryNucleotide)
	return binaryNucleotide


def GroupByLength(bindObj):
	# Group the sequences in the pool by nucleotide length
	nNucs = 0
	maxSeqLen = 0
	bindObj.nSequencesByLengths = np.zeros(bindObj.gene_length, dtype = int)
	with open(("./../data/input/%s" % bindObj.poolfname), 'r') as pool_fp:
		for (l, line) in enumerate(pool_fp):
			if ((l % bindObj.file_skip) == 1):
				nNucs += 1
				seqLen = len(line.strip("\n").strip(" "))
				bindObj.nSequencesByLengths[seqLen] += 1
				if seqLen > maxSeqLen:
					maxSeqLen = seqLen
	bindObj.poolsize = nNucs
	bindObj.max_nuc_length = maxSeqLen
	n_pool_lengths = np.count_nonzero(bindObj.nSequencesByLengths >= 0)
	return n_pool_lengths



def ComputeDerivedParameters(bindObj):
	# Compute parameters required for pattern matching, from the sRNA pool and pstvd sequences.
	# Based on the computed parameters, intialize the respective arrays
	# 1. maximum length of the pstvd sequence
	# 2. total number of nucleotides in the pool.

	# Read the pstvd sequence and compute its length
	with open(("./../data/input/%s" % bindObj.genefname), 'r') as pstvdFid:
		bindObj.gene = pstvdFid.readline().strip("\n")

	bindObj.gene_length = len(bindObj.gene)

	# Group the sequences in the pool by nucleotide length
	n_pool_lengths = GroupByLength(bindObj)

	# Construct the relevant arrays
	bindObj.gene_matrix = np.zeros((4 * bindObj.max_nuc_length, bindObj.gene_length), dtype = int)
	for i in range(bindObj.gene_length):
		bindObj.gene_matrix[:, i] = NucleotideEncoder(StringSlice(bindObj.gene, i, bindObj.max_nuc_length, bindObj.is_circular), bindObj.max_nuc_length, 0)

	bindObj.rawForwMatchFname = ("./../data/output/raw_forward_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	with open(bindObj.rawForwMatchFname, 'w') as rawFid:
		rawFid.write("%d %d\n" % (bindObj.poolsize, bindObj.gene_length))

	bindObj.rawRevMatchFname = ("./../data/output/raw_reverse_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	with open(bindObj.rawForwMatchFname, 'w') as rawFid:
		rawFid.write("%d %d\n" % (bindObj.poolsize, bindObj.gene_length))

	bindObj.forwardMatchFname = ("./../data/output/forward_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	bindObj.forwardMatches = np.zeros((n_pool_lengths, bindObj.gene_length + 1), dtype = np.int32)
	bindObj.reverseMatchFname = ("./../data/output/reverse_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	bindObj.reverseMatches = np.zeros((n_pool_lengths, bindObj.gene_length + 1), dtype = np.int32)

	# Ordering on the non-zero nucleotide lengths
	bindObj.ordering = -1 * np.ones(bindObj.gene_length, dtype = int)
	nzlen = 0
	for i in range(bindObj.gene_length):
		if (bindObj.nSequencesByLengths[i] > 0):
			bindObj.ordering[i] = nzlen
			bindObj.forwardMatches[nzlen, 0] = i
			bindObj.reverseMatches[nzlen, 0] = i
			nzlen = nzlen + 1

	# print("Initialization\nForward\n{}\nReverse\n{}\nOrdering\n{}".format(bindObj.forwardMatches, bindObj.reverseMatches, bindObj.ordering))

	bindObj.forwardMatchCounts = np.zeros(bindObj.gene_length, dtype = np.int32)
	bindObj.reverseMatchCounts = np.zeros(bindObj.gene_length, dtype = np.int32)
	return None


def StringSlice(sequence, start, width, iscirc):
	# Return the substring of a string that starts at a given index and has a given number of characters. The output can be a substring when the string is cyclic.
	substring = []
	for s in range(start, start + width):
		if (iscirc == 0):
			if (s == len(sequence)):
				break
		substring.append(sequence[(s % len(sequence))])
	return substring


def SetBreakPoints(bindObj):
	# Set points which will separate the nucleotides in the sRNA pool into chunks.
	# Each chunk will have at most MAXCONCURRENT nucleotides can be read by a core of the machine.
	nSeqsPartialPool = 0
	# nSeqsRead = 0
	ninvalid = 0
	breakpoints = [[1, 0]]
	with open(("./../data/input/%s" % bindObj.poolfname), 'r') as pool_fp:
		for (l, line) in enumerate(pool_fp):
			sRNASeq = line.strip("\n")
			# every 4th line contains a nucleotide sequence. Others contain data irrelevant to this algorithm
			if (l % bindObj.file_skip == 1):
				if ('N' in sRNASeq):
					ninvalid = ninvalid + 1
				else:
					if (nSeqsPartialPool < bindObj.MAXCONCURRENT):
						nSeqsPartialPool = nSeqsPartialPool + 1
					else:
						breakpoints[-1][-1] = nSeqsPartialPool
						# Start a new pool
						nSeqsPartialPool = 1
						breakpoints.append([l, 0])

		if (breakpoints[-1][-1] == 0):
			breakpoints[-1][-1] = nSeqsPartialPool
	if (ninvalid > 0):
		print("\033[2mThere were %d invalid sequences -- having \"N\" in them.\033[0m" % (ninvalid))
	return breakpoints


def PartialBinding(pool_matrix, gene_matrix, gene_length, length_offsets, pool_lengths, ordering, matches, counts):
	# Compute the frequency of matchings in the sRNA pool with every continuous PSTVd subsequence.
	# First, compute the Pattern matchings by matrix multiplication of the sRNA sequence encoding matrix with the PSTVd sequences encoding matrix.
	# If there was a perfect matching between a PSTVd subsequene and the sRNA nucleotide in the pool, the corresponding element generated in the product will be equal to the length of the sRNA nucleotide.
	# The amount by which the product element is lesser than the length of the sRNA nucleotide, is precisely the number of mismatches.
	# If the product doescribes matching upto the right tolerence, then reduce the product to 1, or else to 0.
	# matchings = matrix whose i,j entry is equal to the number of locations at which the ith nucleotide in the sRNA pool matched with the pstvd sequence of length L starting at index j.
	matchings = np.dot(pool_matrix, gene_matrix)
	# reduced = matrix whose i,j entry is 1 if the ith nucleotide matched with the pstvd sequence (of length L) starting at index j.
	reduced = np.greater_equal(matchings, length_offsets)
	(pool_indices, gene_indices) = np.nonzero(reduced)
	for i in range(pool_indices.shape[0]):
		matches[ordering[pool_lengths[pool_indices[i]]] * gene_length + gene_indices[i]] += 1
		counts[pool_lengths[pool_indices[i]]] += 1
	return None


def	PreparePartialBinding(n_seqs, skip_lines, reverse, rnap):
	# Partially load the sRNA pool with either the remaining or the next MAXCONCURRENT nucleotides, from the main pool, whichever is lesser.
	# The number of nucleotide sequences to be read can be less than the above estimate -- because valid nucleotide sequences should not have 'N' in them.
	pool_matrix = np.zeros((n_seqs, 4 * rnap.max_nuc_length), dtype = int)
	length_offsets = np.zeros((n_seqs, rnap.gene_length), dtype = int)
	pool_lengths = np.zeros(n_seqs, dtype = int)
	
	with open(("./../data/input/%s" % rnap.poolfname), 'r') as pool_fp:
		nSkipped = 0
		nRead = 0
		completed = 1
		
		# Skip the nucleotides that have been read previously, i.e, those which appear before skip_lines.
		for line_number in range(skip_lines):
			pool_fp.readline()

		# print("Skipped {} lines, line_number: {}.".format(skip_lines, line_number))

		while (nRead < n_seqs):
			line_number += 1
			line = pool_fp.readline()
			# Every k-th (k = file_skip) line contains a nucleotide sequence.
			# Others contain data irrelevant to this algorithm.
			if ((line_number % rnap.file_skip) == 1):
				sRNASeq = line.strip("\n").strip(" ")
				if ('N' in sRNASeq):
					pass
				else:
					pool_matrix[nRead, :] = NucleotideEncoder(sRNASeq, rnap.max_nuc_length, reverse)
					pool_lengths[nRead] = len(sRNASeq)
					length_offsets[nRead, :] = (len(sRNASeq) - rnap.tolerance) * np.ones(rnap.gene_length, dtype = int)
					nRead += 1
	return (pool_matrix, length_offsets, pool_lengths)


def RecordBinding(rnap):
	# Save the matching frequency array
	WriteMatrixToFile(rnap.forwardMatchFname, rnap.forwardMatches, 0, 0, 0)
	# print("Reverse matches\n{}".format(rnap.reverseMatches))
	WriteMatrixToFile(rnap.reverseMatchFname, rnap.reverseMatches, 0, 0, 0)
	print("\033[2mResults written to %s and %s.\033[0m" % (rnap.forwardMatchFname, rnap.reverseMatchFname))
	return None


def LogBinding(rnap):
	# Store the details of the pattern matching run into a log file. The file is to be appened.
	topology = ["linear", "circular"]
	with open("log.txt", 'a') as logfp:
		logfp.write("\n***********************************************************************\n")
		logfp.write("Details for simulations done on %s at %s.\n" % tuple(rnap.dateStamp.split(" ")))

		logfp.write("Length of the PSTVd sequence in %s: %d.\n" % (rnap.genefname, rnap.gene_length))
		logfp.write("Number of nucleotides in the sRNA pool in %s: %d.\n" % (rnap.poolfname, rnap.poolsize))
		logfp.write("Number of mismatches allowed: %d.\n" % rnap.tolerance)
		logfp.write("Topology of matching: %s.\n" % topology[rnap.is_circular])
		
		logfp.write("================================================================\n")

		logfp.write("{:^6} | {:^8} | {:^8} | {:^8} | {:^8} | {:^8}\n".format("Length", "In pool", "Forward", "Reverse", "Total", "% mapped"))
		logfp.write("----------------------------------------------------------------\n")
		total = 0
		for l in range(rnap.gene_length):
			if (rnap.nSequencesByLengths[l] > 0):
				total += rnap.nSequencesByLengths[l]
				total_matches = rnap.forwardMatchCounts[l] + rnap.reverseMatchCounts[l]
				percentage_mapped = total_matches/rnap.nSequencesByLengths[l] * 100
				logfp.write("{:^6} | {:^8} | {:^8} | {:^8} | {:^8} | {:^8}\n".format(l, rnap.nSequencesByLengths[l], rnap.forwardMatchCounts[l], rnap.reverseMatchCounts[l], total_matches, "%.2f" % (percentage_mapped)))
		logfp.write("----------------------------------------------------------------\n")
		total_matches = np.sum(rnap.forwardMatchCounts) + np.sum(rnap.reverseMatchCounts)
		percentage_mapped = total_matches/rnap.poolsize * 100
		logfp.write("{:^6} | {:^8} | {:^8} | {:^8} | {:^8} | {:^8}\n".format("", total, np.sum(rnap.forwardMatchCounts), np.sum(rnap.reverseMatchCounts), total_matches, "%.2f" % percentage_mapped))
		logfp.write("================================================================\n")

		logfp.write("Total simulation time: %d seconds." % rnap.runtime)
		logfp.write("\n***********************************************************************\n")
	return None


def WriteMatrixToFile(fname, mat, appendMode, sparse = 0, binary = 0):
	# Write a matrix to file
	(nRows, nCols) = mat.shape
	if appendMode == 0:
		if os.path.isfile(fname):
			os.remove(fname)
	with open(fname, 'a') as matFid:
		if appendMode == 0:
			matFid.write("%d %d\n" % (nRows, nCols))
		for ri in range(nRows):
			for ci in range(nCols):
				if ((sparse == 1) and (binary == 1)):
					if mat[ri, ci] == 0:
						pass
					else:
						matFid.write("%d " % ci)
				else:
					if (sparse == 1):
						if mat[ri, ci] == 0:
							pass
						else:
							matFid.write("%d %d" % (ci, mat[ri, ci]))
					else:
						matFid.write("%d " % (mat[ri, ci]))
			matFid.write("\n")
	return None


def ShowOutput(rnap, ncores):
	# Print the output of the matchings to the console.
	print("_/ Matching was completed in %d seconds on each of %d cores." % (rnap.runtime, ncores))
	# print("Length\tIn pool\tFoward\tReverse\tTotal")
	print("=====================================================================")
	print("{:^6} | {:^8} | {:^8} | {:^8} | {:^8} | {:^8}".format("Length", "In pool", "Forward", "Reverse", "Total", "Reads per 100"))
	print("---------------------------------------------------------------------")
	total = 0
	for l in range(rnap.gene_length):
		if (rnap.nSequencesByLengths[l] > 0):
			# print("%d\t%d\t%d\t%d\t%d" % (l, rnap.nSequencesByLengths[l], rnap.forwardMatchCounts[l], rnap.reverseMatchCounts[l], rnap.forwardMatchCounts[l] + rnap.reverseMatchCounts[l]))
			total += rnap.nSequencesByLengths[l]
			total_matches = rnap.forwardMatchCounts[l] + rnap.reverseMatchCounts[l]
			percentage_mapped = total_matches/rnap.nSequencesByLengths[l] * 100
			print("{:^6} | {:^8} | {:^8} | {:^8} | {:^8} | {:^8}".format(l, rnap.nSequencesByLengths[l], rnap.forwardMatchCounts[l], rnap.reverseMatchCounts[l], total_matches, "%.2f" % (percentage_mapped)))
	print("----------------------------------------------------------------")
	total_matches = np.sum(rnap.forwardMatchCounts) + np.sum(rnap.reverseMatchCounts)
	percentage_mapped = total_matches/rnap.poolsize * 100
	print("{:^6} | {:^8} | {:^8} | {:^8} | {:^8} | {:^8}".format("", total, np.sum(rnap.forwardMatchCounts), np.sum(rnap.reverseMatchCounts), total_matches, "%.2f" % percentage_mapped))
	print("=====================================================================")
	return None


def RunMatchings(rnap, ncores):
	# Run the matchings with the inputs loaded into the object
	topology = ["linear", "circular"]
	print("Now Running\n\tgene: %s\n\tpool: %s\n\ttopology: %s\n\ttolerance: %d\n\tcores: %d" % (rnap.genefname, rnap.poolfname, topology[rnap.is_circular], rnap.tolerance, ncores))
	
	breakpoints = SetBreakPoints(rnap)
	nBreaks = len(breakpoints)

	start = time.time()
	for d in range(2):
		# change above 'range (2)' to have both forward and reverse matching
		if (d == 0):
			print("\033[2mForward matchings:")
		else:
			print("\033[2mReverse matchings:")
		
		# Define the counts and the matches arrays
		matches = Array('i', rnap.forwardMatches.size - rnap.forwardMatches.shape[0])
		counts = Array('i', rnap.gene_length)
		
		completed = 0

		while (completed < nBreaks):
			local_start = time.time()
			launch = min(nBreaks - completed, ncores)

			# Reset the counts and matches array.
			for i in range(rnap.forwardMatchCounts.size):
				counts[i] = 0
			for i in range(rnap.forwardMatches.size - rnap.forwardMatches.shape[0]):
				matches[i] = 0

			processes = []
			if ((QUIET == 0)):
				if (launch == 0):
					sys.stderr.write("ncores = 0\n")

			for p in range(launch):
				lines_to_skip = breakpoints[completed + p][0]
				seqs_to_read = breakpoints[completed + p][1]
				(pool_matrix, length_offsets, pool_lengths) = PreparePartialBinding(seqs_to_read, lines_to_skip, d, rnap)
				processes.append(Process(target = PartialBinding, args = (pool_matrix, rnap.gene_matrix, rnap.gene_length, length_offsets, pool_lengths, rnap.ordering, matches, counts)))
				# PartialBinding(pool_matrix, rnap.gene_matrix, rnap.gene_length, length_offsets, pool_lengths, rnap.ordering, matches, counts)

			for p in range(launch):
				processes[p].start()

			for p in range(launch):
				processes[p].join()

			# Gathering results
			if (d == 0):
				for i in range(rnap.forwardMatches.shape[0]):
					for j in range(rnap.gene_length):
						rnap.forwardMatches[i, j + 1] = rnap.forwardMatches[i, j + 1] + matches[i * rnap.gene_length + j]
				for j in range(rnap.gene_length):
					rnap.forwardMatchCounts[j] = rnap.forwardMatchCounts[j] + counts[j]
			else:
				for i in range(rnap.reverseMatches.shape[0]):
					for j in range(rnap.gene_length):
						rnap.reverseMatches[i, j + 1] = rnap.reverseMatches[i, j + 1] + matches[i * rnap.gene_length + j]
				for j in range(rnap.gene_length):
					rnap.reverseMatchCounts[j] = rnap.reverseMatchCounts[j] + counts[j]

			completed = completed + launch

			local_runtime = time.time() - local_start
			# Predict the remaining time to finish all the matching: remaining = time taken for single iteration x number of iterations left
			# The time taken for a single iteration can be measured in the first iterationa and the number of iterations left is at most the total number of sequences in the pool divided by the number of sequences read in a chunk.
			remaining = local_runtime * (nBreaks - completed)

			if (QUIET == 0):
				print("\033[2m ......... %d%% done, about %d minutes %d seconds left\033[0m" % (100 * completed/float(nBreaks), (remaining/60), (int(remaining) % 60)))
			else:
				print("\r\033[2m ......... %d%% done, about %d minutes %d seconds left\033[0m" % (100 * completed/float(nBreaks), (remaining/60), (int(remaining) % 60))),
				sys.stdout.flush()

	rnap.runtime = (time.time() - start)
	rnap.dateStamp = time.strftime("%d/%m/%Y %H:%M:%S")
	return None


def Bind(rnap):
	# Solve the Viroid binding problem.
	# This is a wrapper function for the Python package.
	# The input is an object of the RNAProfiler class. The attributes:
	# genefname, poolfname, is_circular and tolerance should be loaded before this function call.
	
	# Compute derived parameters
	ComputeDerivedParameters(rnap)
	
	# Run pattern matching
	RunMatchings(rnap, ncores)

	# Save the matching information to a file
	RecordBinding(rnap)
	LogBinding(rnap)
	return None


if __name__ == '__main__':
	
	# Display the logo and license information
	DisplayLogoLicense()
	# Check if all the required packages exist
	CheckDependencies()

	# Complete the binding mechanism by pattern matching batches of the sRNA nucleotides with the PSTVd sequence. The batches are taken such that no batch is greater in size than MAXCONCURRENT and all the batches together constitute all the nucleotides in the sRNA pool.
	if (len(sys.argv) < 2):
		Usage()
		exit(0)
	else:
		if (len(sys.argv) > 2):
			task = int(sys.argv[2].strip("\n"))
		else:
			task = -1
		input_file = ("./../data/input/%s" % sys.argv[1].strip("\n"))


	# Reading the gene, pool files and tolerance.
	if not os.path.isfile(input_file):
		print("\033[2mError: Input file %s does not exist in data/input.\033[0m" % (input_file))
		exit(0);
	else:
		input_fp = open(input_file, "r")

	n_task = 0
	stop = 0

	while (stop == 0):
		line = input_fp.readline().strip("\n").strip(" ")
		if not (line[0] == "#"):
			if ("quit" in line):
				stop = 1
				continue
			n_task += 1
			if ((n_task == task) or (task == -1)):
				linecontents = list(map(lambda ln: ln.strip("\n").strip(" "), line.split(" ")))
				# Initialization of the RNAProfiler object and its attributes
				rnap = RNAProfiler()
				rnap.genefname = linecontents[0]
				rnap.poolfname = linecontents[1]
				rnap.file_skip = int(linecontents[2])
				rnap.tolerance = int(linecontents[3])
				rnap.is_circular = int(linecontents[4])
				ncores = 1
				# ncores = int(linecontents[5])
				
				# If the task is not specified, do all the tasks in the input file.
				if (task > -1):
					stop = 1

				# Compute derived parameters
				ComputeDerivedParameters(rnap)
	
				# Run pattern matching
				RunMatchings(rnap, ncores)

				# Print output
				ShowOutput(rnap, ncores)

				# Save the matching information to a file
				RecordBinding(rnap)
				LogBinding(rnap)

	input_fp.close()