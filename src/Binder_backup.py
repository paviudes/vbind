import sys
from multiprocessing import Process, Array
import os
import time
import numpy as np

DEBUG = 0
QUIET = 0

class Binder():
	"""
	Binder class: Contains all the information required to perform a round of Pattern matching
	string types
	poolfname = name of the file containing the nucleotide sequences
	genefname = name of the file containing the pstvd sequences
	gene = pstvd sequence
	dateStamp = date and time when the pattern matching simulation was completed
	rawForwMatchFname = name of the file containing the raw forward matching data -- matrix (i,j) entry is 1 iff ith sRNA sequence matched with in the forward alignment with the PSTVd subsequence starting at index j.
	rawRevMatchFname =  name of the file containing the raw reverse matching data -- matrix (i,j) entry is 1 iff ith sRNA sequence matched with in the reverse alignment with the PSTVd subsequence starting at index j.
	forwardMatchFname = name of the file containing the matching data -- matrix whose (i,j) entry is the number of times a sequence of length j forward-matched with a PSTVd subsequence starting at index j.
	reverseMatchFname = name of the file containing the matching data -- matrix whose (i,j) entry is the number of times a sequence of length j reverse-matched with a PSTVd subsequence starting at index j.
	forwMappingFname = name of the file containing the forward mapping data
	revMappingFname = name of the file containing the reverse mapping data
	integer types
	max_nuc_length = maximum length of a nucleotide sequence
	gene_length = length of the pstvd sequence
	MAXCONCURRENT = 100,000
	poolsize = number of nucleotides (of all lengths) in the sRNA pool. The sRNA pool is the one that is currently being matched for. This will be limited to a maximum of MAXCONCURRENT. To handle a pool of a larger size, we must repeat the matching separately for smaller pools and gather the matching data.
	tolerance = maximum number of mismatches authorized to validate a matching.
	floating types
	runtime = total running time of the pattern matching algorithm
	Integer arrays
	breakpoints = nucleotide indices and number of nucleotides that define a chunk
	gene_matrix = matrix whose column j is the binary encoding of the pstvd sequence of length L starting at position j. L is the maximum length of a nucleotide in the current sRNA pool.
	sRNAPoolMatrix = matrix whose jth row is the binary encoding of the jth nucleotide in the sRNA pool.
	sRNASeqLengths = matrix whose ith row contains the length of the ith sRNA nucleotide (minus the tolerance) repeated as many times as the length of the length of the PSTVd sequence
	nSequencesByLengths = number of nucleotide sequences in the pool of every length.
	forwardMatches = matrix whose i,j entry denotes the number of nucleotides of length j that matched in the forward alignment direction with the pstvd substring starting at index i.
	reverseMatches = matrix whose i,j entry denotes the number of nucleotides of length j that matched in the reverse alignment direction with the pstvd substring starting at index i.
	forwardMapping = vector whose ith index specifies how many times the ith character in the PSTVd sequence appeared in a forward matching with any nucleotide
	reverseMapping = vector whose ith index specifies how many times the ith character in the PSTVd sequence appeared in a reverse matching with any nucleotide
	"""

	def __init__(self):
		# Initialize all the variables in the scope of the Binder class
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
		self.isCircular = 1

		self.gene_matrix = np.array([])
		self.forwardMatches = np.array([])
		self.reverseMatches = np.array([])
		self.forwardMatchCounts = np.array([])
		self.reverseMatchCounts = np.array([])
		return None


# Load a Binder object
def Load(bindObj, interfaceFile = "NAN"):
	# Load all the information required in the Binder object for Pattern matching
	# Read the interface file and load the parameters for Binding
	if os.path.isfile(interfaceFile):
		bindObj.interfaceFile = interfaceFile
		if (QUIET == 0):
			print("\033[94mLoading information required for pattern matching from %s.\033[0m" % bindObj.interfaceFile)

		with open(bindObj.interfaceFile, 'r') as uiFid:
			# 1. file name of pstvd sequence
			# 2. file name of nucleotides in the pool
			# 3. tolerance
			bindObj.genefname = uiFid.readline().strip("\n").strip(" ")
			bindObj.poolfname = uiFid.readline().strip("\n").strip(" ")
			bindObj.tolerance = int(uiFid.readline().strip("\n").strip(" "))
	else:
		LoadBinderFromConsole(bindObj)
	return None


def LoadBinderFromConsole(bindObj):
	# Load a Binder object from taking console inputs
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
		print("\033[95mNucleotide sequence: %s and encoding produced for direction %d:" % (nucleotide, reverse))
		print(binaryNucleotide)
	return binaryNucleotide


def ComputeDerivedParameters(bindObj):
	# Compute parameters required for pattern matching, from the sRNA pool and pstvd sequences.
	# Based on the computed parameters, intialize the respective arrays
	# 1. maximum length of the pstvd sequence
	# 2. total number of nucleotides in the pool.

	# Read the pstvd sequence and compute its length
	with open(("./../data/input/%s" % bindObj.genefname), 'r') as pstvdFid:
		bindObj.gene = pstvdFid.readline().strip("\n")

	bindObj.gene_length = len(bindObj.gene)

	# Compute the total number of nucleotides and the maximum length of a nucleotide sequence
	nNucs = 0
	maxSeqLen = 0
	bindObj.nSequencesByLengths = np.zeros(bindObj.gene_length, dtype = int)
	with open(("./../data/input/%s" % bindObj.poolfname), 'r') as sRNAFid:
		for (lNo, line) in enumerate(sRNAFid):
			if ((lNo % 2) == 1):
				print("line = {}".format(line))
				nNucs = nNucs + 1
				seqLen = len(line.strip("\n").strip(" "))
				bindObj.nSequencesByLengths[seqLen] = bindObj.nSequencesByLengths[seqLen] + 1
				if seqLen > maxSeqLen:
					maxSeqLen = seqLen
	bindObj.poolsize = nNucs
	bindObj.max_nuc_length = maxSeqLen
	print("Sequence read")

	# Construct the relevant arrays
	bindObj.gene_matrix = np.zeros((4 * bindObj.max_nuc_length, bindObj.gene_length), dtype = int)
	for i in range(bindObj.gene_length):
		bindObj.gene_matrix[:, i] = NucleotideEncoder(StringSlice(bindObj.gene, i, bindObj.max_nuc_length, bindObj.isCircular), bindObj.max_nuc_length, 0)

	bindObj.rawForwMatchFname = ("./../data/output/raw_forward_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	with open(bindObj.rawForwMatchFname, 'w') as rawFid:
		rawFid.write("%d %d\n" % (bindObj.poolsize, bindObj.gene_length))

	bindObj.rawRevMatchFname = ("./../data/output/raw_reverse_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	with open(bindObj.rawForwMatchFname, 'w') as rawFid:
		rawFid.write("%d %d\n" % (bindObj.poolsize, bindObj.gene_length))

	bindObj.forwardMatchFname = ("./../data/output/forward_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	bindObj.forwardMatches = np.zeros((np.count_nonzero(bindObj.nSequencesByLengths), bindObj.gene_length + 1), dtype = int)
	bindObj.reverseMatchFname = ("./../data/output/reverse_matchings_%s_%s_tol%d.txt" % (bindObj.poolfname[:-4], bindObj.genefname[:-4], bindObj.tolerance))
	bindObj.reverseMatches = np.zeros((np.count_nonzero(bindObj.nSequencesByLengths), bindObj.gene_length + 1), dtype = int)

	# Ordering on the non-zero nucleotide lengths
	bindObj.ordering = (-1) * np.ones(bindObj.gene_length, dtype = int)
	nzlen = 0
	for i in range(bindObj.gene_length):
		if (bindObj.nSequencesByLengths[i] > 0):
			bindObj.ordering[i] = nzlen
			bindObj.forwardMatches[nzlen, 0] = i
			bindObj.reverseMatches[nzlen, 0] = i
			nzlen = nzlen + 1

	# print("Initialization\nForward\n{}\nReverse\n{}".format(bindObj.forwardMatches, bindObj.reverseMatches))

	bindObj.forwardMatchCounts = np.zeros(bindObj.gene_length, dtype = int)
	bindObj.reverseMatchCounts = np.zeros(bindObj.gene_length, dtype = int)
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
	nSeqsRead = 0
	ninvalid = 0
	breakpoints = [[0, 0, 0]]
	with open(("./../data/input/%s" % bindObj.poolfname), 'r') as sRNAFid:
		for (lNo, line) in enumerate(sRNAFid):
			sRNASeq = line.strip("\n")
			# every 4th line contains a nucleotide sequence. Others contain data irrelevant to this algorithm
			if (lNo % 4 == 1):
				if ('N' in sRNASeq):
					ninvalid = ninvalid + 1
				else:
					if (nSeqsPartialPool < bindObj.MAXCONCURRENT):
						nSeqsPartialPool = nSeqsPartialPool + 1
						nSeqsRead = nSeqsRead + 1
					else:
						breakpoints[-1][-1] = nSeqsPartialPool
						# Start a new pool
						nSeqsPartialPool = 0
						breakpoints.append([lNo, nSeqsRead, 0])
		if (breakpoints[-1][-1] == 0):
			breakpoints[-1][-1] = nSeqsPartialPool
	# print("nSeqsPartialPool = {}, nSeqsRead = {}, ninvalid = {}, breakpoints = {}".format(nSeqsPartialPool, nSeqsRead, ninvalid, breakpoints))
	return breakpoints


def PartialBinding(bindObj, nSkip, startSeq, numSeqs, reverse, matches, counts):
	# Partially load the sRNA pool with either the remaining or the next MAXCONCURRENT nucleotides, from the main pool, whichever is lesser.
	# The number of nucleotide sequences to be read can be less than the above estimate -- because valid nucleotide sequences should not have 'N' in them.
	sRNAPoolMatrix = np.zeros((numSeqs, 4 * bindObj.max_nuc_length), dtype = int)
	seqLenOffsets = np.zeros((numSeqs, bindObj.gene_length), dtype = int)
	sRNASeqLengths = np.zeros(numSeqs, dtype = int)
	directions = ["forward", "reverse"]

	with open(("./../data/input/%s" % bindObj.poolfname), 'r') as sRNAFid:
		nSkipped = 0
		nRead = 0
		completed = 1
		for (lNo, line) in enumerate(sRNAFid):
			# every 4th line contains a nucleotide sequence. Others contain data irrelevant to this algorithm
			# Skip the nucleotides that have been read previously, i.e, those which appear before seqStartLineNo
			# We need to skip seqStart*4 lines.
			if (lNo < nSkip):
				continue
			else:
				if ((lNo % 4) == 1):
					sRNASeq = line.strip("\n").strip(" ")
					# print("Read\n{}".format(sRNASeq))
					if ('N' in sRNASeq):
						pass
					else:
						if (nRead < numSeqs):
							sRNAPoolMatrix[nRead, :] = NucleotideEncoder(sRNASeq, bindObj.max_nuc_length, reverse)
							sRNASeqLengths[nRead] = len(sRNASeq)
							seqLenOffsets[nRead, :] = (len(sRNASeq) - bindObj.tolerance) * np.ones(bindObj.gene_length, dtype = int)
							nRead = nRead + 1

	# Compute the frequency of matchings in the sRNA pool with every continuous PSTVd subsequence.
	# First, compute the Pattern matchings by matrix multiplication of the sRNA sequence encoding matrix with the PSTVd sequences encoding matrix.
	# If there was a perfect matching between a PSTVd subsequene and the sRNA nucleotide in the pool, the corresponding element generated in the product will be equal to the length of the sRNA nucleotide.
	# The amount by which the product element is lesser than the length of the sRNA nucleotide, is precisely the number of mismatches.
	# If the product doescribes matching upto the right tolerence, then reduce the product to 1, or else to 0.

	# matchings = matrix whose i,j entry is equal to the number of locations at which the ith nucleotide in the sRNA pool matched with the pstvd sequence of length L starting at index j.
	matchings = np.dot(sRNAPoolMatrix, bindObj.gene_matrix)
	# reducedToTolerance = matrix whose i,j entry is 1 if the ith nucleotide matched with the pstvd sequence (of length L) starting at index j.
	reducedToTolerance = np.greater_equal(matchings, seqLenOffsets)
	(localNucIndices, pstvdIndices) = np.nonzero(reducedToTolerance)
	for i in range(localNucIndices.shape[0]):
		matches[bindObj.ordering[sRNASeqLengths[localNucIndices[i]]] * bindObj.gene_length + pstvdIndices[i]] += 1
		counts[sRNASeqLengths[localNucIndices[i]]] += 1

	# print("sRNAPoolMatrix\n{}".format(sRNAPoolMatrix))
	# print("bindObj.gene_matrix\n{}".format(bindObj.gene_matrix))
	# print("matchings\n{}".format(matchings))
	# print("reducedToTolerance\n{}".format(np.where(reducedToTolerance, 1, 0)))
	return None


def RecordBinding(bindObj):
	# Save the matching frequency array
	WriteMatrixToFile(bindObj.forwardMatchFname, bindObj.forwardMatches, 0, 0, 0)
	# print("Reverse matches\n{}".format(bindObj.reverseMatches))
	WriteMatrixToFile(bindObj.reverseMatchFname, bindObj.reverseMatches, 0, 0, 0)
	return None


def LogBinding(bindObj):
	# Store the details of the pattern matching run into a log file. The file is to be appened.
	logFile = "details.txt"
	with open(logFile, 'w') as logFid:
		logFid.write("\n***************************************\n")

		logFid.write("Log file for simulations done on, at, %s.\n" % bindObj.dateStamp)

		logFid.write("Length of the PSTVd sequence in %s: %d.\n" % (bindObj.genefname, bindObj.gene_length))
		logFid.write("Number of nucleotides in the sRNA pool in %s: %d.\n" % (bindObj.poolfname, bindObj.poolsize))
		logFid.write("Number of mismatches allowed: %d.\n" % bindObj.tolerance)
		logFid.write("Length\tIn pool\tFoward\tReverse\n")
		for li in range(bindObj.gene_length):
			if (bindObj.nSequencesByLengths[li] > 0):
				logFid.write("%d\t%d\t%d\t%d\n" % (li, bindObj.nSequencesByLengths[li], bindObj.forwardMatchCounts[li], bindObj.reverseMatchCounts[li]))
		logFid.write("-------\n")
		logFid.write("Total:\t%d\t%d\t%d\n\n" % (bindObj.poolsize, np.sum(bindObj.forwardMatchCounts), np.sum(bindObj.reverseMatchCounts)))

	with open(logFile, 'a') as logFid:
		logFid.write("\n\n")
		logFid.write("Total simulation time: %d seconds." % bindObj.runtime)
		logFid.write("\n***************************************\n")
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


if __name__ == '__main__':
	# Complete the binding mechanism by pattern matching batches of the sRNA nucleotides with the PSTVd sequence. The batches are taken such that no batch is greater in size than MAXCONCURRENT and all the batches together constitute all the nucleotides in the sRNA pool.
	# Initialization of the Binder object and its attributes
	#print("Received signal: {}".format(sys.argv))
	rnas = Binder()
	# Reading the gene, pool files and tolerance.
	if os.path.isfile(sys.argv[1], "r"):
		with open(sys.argv[1], "r") as fp:
			for (l, line) in enumerate(fp):
				if (l == int(sys.argv[2])):
					# print("Found relevant input in line: {}".format(line))
					linecontents = list(map(lambda ln: ln.strip("\n").strip(" "), line.split(" ")))
					# print("linecontents = {}".format(linecontents))
					rnas.genefname = linecontents[0]
					rnas.poolfname = linecontents[1]
					rnas.tolerance = int(linecontents[2])
					nCores = int(linecontents[3])
	ComputeDerivedParameters(rnas)
	breakpoints = SetBreakPoints(rnas)
	nBreaks = len(breakpoints)

	print("{} breakpoints\n{}".format(nBreaks, breakpoints))

	print("\033[93mMatching with %d sRNA nucleotides and a gene of length %d. Maximum length of a nucleotide sequence is %d.\033[0m" % (rnas.poolsize, rnas.gene_length, rnas.max_nuc_length))
	startTime = time.time()
	matches = Array('i', rnas.forwardMatches.size - rnas.forwardMatches.shape[0])
	counts = Array('i', rnas.forwardMatchCounts.size)
	for di in range(2):
		# change above 'range (2)' to have both forward and reverse matching
		if (di == 0):
			print("\033[94mForward matchings:")
		else:
			print("\033[94mReverse matchings:")
		completed = 0
		while (completed < nBreaks):
			start = time.time()
			launchNow = min(nBreaks - completed, nCores)

			# Reset the counts and the matches arrays
			for i in range(rnas.forwardMatchCounts.size):
				counts[i] = 0
			for i in range(rnas.forwardMatches.size - rnas.forwardMatches.shape[0]):
				matches[i] = 0

			processes = []
			if ((QUIET == 0)):
				if (launchNow == 0):
					sys.stderr.write("nCores = 0\n")

			for pi in range(launchNow):
				linesToSkip = breakpoints[completed + pi][0]
				seqsReadSoFar = breakpoints[completed + pi][1]
				seqsToRead = breakpoints[completed + pi][2]
				processes.append(Process(target = PartialBinding, args = (rnas, linesToSkip, seqsReadSoFar, seqsToRead, di, matches, counts)))

			for pi in range(launchNow):
				processes[pi].start()

			for pi in range(launchNow):
				processes[pi].join()

			# Gathering results
			if (di == 0):
				# for i in range(rnas.forwardMatches.shape[0]):
				#       for j in range(rnas.gene_length):
				#               print("matches[%d][%d] = %d" % (i, j, matches[i * (rnas.forwardMatches.shape[1] - 1) + j]))
				for i in range(rnas.forwardMatches.shape[0]):
					for j in range(rnas.gene_length):
						rnas.forwardMatches[i, j + 1] += matches[i * (rnas.forwardMatches.shape[1] - 1) + j]
				for i in range(rnas.gene_length):
					rnas.forwardMatchCounts[i] += counts[i]
			else:
				for i in range(rnas.reverseMatches.shape[0]):
					for j in range(rnas.gene_length):
						rnas.reverseMatches[i, j + 1] += matches[i * (rnas.reverseMatches.shape[1] - 1) + j]
				for i in range(rnas.gene_length):
					rnas.reverseMatchCounts[i] += counts[i]

			completed = completed + launchNow

			iterationEndTime = time.time()
			runtime = iterationEndTime - start
			# Predict the remaining time to finish all the matching: remaining = time taken for single iteration x number of iterations left
			# The time taken for a single iteration can be measured in the first iterationa and the number of iterations left is at most the total number of sequences in the pool divided by the number of sequences read in a chunk.
			remaining = runtime * (nBreaks - completed)

			if (QUIET == 0):
				print("\033[95m ......... %d%% done, about %d minutes %d seconds left\033[0m" % (100 * completed/float(nBreaks), (remaining/60), (int(remaining) % 60)))
			else:
				print("\r\033[95m ......... %d%% done, about %d minutes %d seconds left\033[0m" % (100 * completed/float(nBreaks), (remaining/60), (int(remaining) % 60))),
				sys.stdout.flush()

		if (QUIET > 0):
			print("")

	endTime = time.time()
	rnas.runtime = (endTime - startTime)
	rnas.dateStamp = time.strftime("%d/%m/%Y %H:%M:%S")

	if (QUIET == 0):
		print("\033[92m_/ Matching was completed in %d seconds on each of %d cores.\033[0m" % (rnas.runtime, nCores))
	else:
		print("\r\033[92m_/ Matching was completed in %d seconds on each of %d cores.\033[0m" % (rnas.runtime, nCores))

	print("\033[92mSummary:\033[0m")
	print("\033[92mLength\tIn pool\tFoward\tReverse\033[0m")
	for li in range(rnas.gene_length):
		if (rnas.nSequencesByLengths[li] > 0):
			print("\033[92m%d\t%d\t%d\t%d\033[0m" % (li, rnas.nSequencesByLengths[li], rnas.forwardMatchCounts[li], rnas.reverseMatchCounts[li]))
	print("\033[92m--------------\033[0m")
	print("\033[92mTotal:\t%d\t%d\t%d\033[0m" % (rnas.poolsize, np.sum(rnas.forwardMatchCounts), np.sum(rnas.reverseMatchCounts)))

	print("Forward matchings\n{}".format(rnas.forwardMatches))
	print("Reverse matchings\n{}".format(rnas.reverseMatches))

	# Save the matching information to a file
	RecordBinding(rnas)
	LogBinding(rnas)

	print("\033[93mSummary written to %s and %s.\033[0m" % (rnas.forwardMatchFname, rnas.reverseMatchFname))