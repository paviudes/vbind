import sys
from multiprocessing import Queue, Process
import os
import time
import numpy as np

DEBUG = 0
QUIET = 0

class Binder():
	"""
	Binder class: Contains all the information required to perform a round of Pattern matching

	string types
	sRNAPoolFname = name of the file containing the nucleotide sequences
	pstvdFname = name of the file containing the pstvd sequences
	pstvdSeq = pstvd sequence
	dateStamp = date and time when the pattern matching simulation was completed
	rawForwMatchFname = name of the file containing the raw forward matching data -- matrix (i,j) entry is 1 iff ith sRNA sequence matched with in the forward alignment with the PSTVd subsequence starting at index j.
	rawRevMatchFname =  name of the file containing the raw reverse matching data -- matrix (i,j) entry is 1 iff ith sRNA sequence matched with in the reverse alignment with the PSTVd subsequence starting at index j.
	forwardMatchFname = name of the file containing the matching data -- matrix whose (i,j) entry is the number of times a sequence of length j forward-matched with a PSTVd subsequence starting at index j.
	reverseMatchFname = name of the file containing the matching data -- matrix whose (i,j) entry is the number of times a sequence of length j reverse-matched with a PSTVd subsequence starting at index j.

	forwMappingFname = name of the file containing the forward mapping data
	revMappingFname = name of the file containing the reverse mapping data

	integer types
	maxNucleotideLength = maximum length of a nucleotide sequence
	pstvdLen = length of the pstvd sequence
	MAXCONCURRENT = 100,000
	totalNucleotides = number of nucleotides (of all lengths) in the sRNA pool. The sRNA pool is the one that is currently being matched for. This will be limited to a maximum of MAXCONCURRENT. To handle a pool of a larger size, we must repeat the matching separately for smaller pools and gather the matching data.
	tolerance = maximum number of mismatches authorized to validate a matching.
	
	floating types
	runtime = total running time of the pattern matching algorithm

	Integer arrays
	breakpoints = nucleotide indices and number of nucleotides that define a chunk
	pstvdMatrix = matrix whose column j is the binary encoding of the pstvd sequence of length L starting at position j. L is the maximum length of a nucleotide in the current sRNA pool.
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
		self.sRNAPoolFname = "NAN"
		self.pstvdFname = "NAN"
		self.matchFrequencyFname = "NAN"
		self.dateStamp = ""
		
		self.rawForwMatchFname = ""
		self.rawRevMatchFname = ""
		
		self.forwardMatchFname = ""
		self.reverseMatchFname = ""

		self.runtime = 0

		self.pstvdLen = 0
		self.MAXCONCURRENT = 100000000
		self.totalNucleotides = 0
		self.tolerance = 0
		self.maxNucleotideLength = 0

		self.pstvdMatrix = np.array([])
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
			bindObj.pstvdFname = uiFid.readline().strip("\n").strip(" ")
			bindObj.sRNAPoolFname = uiFid.readline().strip("\n").strip(" ")
			bindObj.tolerance = int(uiFid.readline().strip("\n").strip(" "))
	else:
		LoadBinderFromConsole(bindObj)
	return None


def LoadBinderFromConsole(bindObj):
	# Load a Binder object from taking console inputs
	print("\033[93mEnter the name of the file containing the sRNA pool.\033[0m")
	print(">>"),
	bindObj.sRNAPoolFname = raw_input().strip("\n")

	print("\033[93mEnter the name of the file containing the PSTVd sequence.\033[0m")
	print(">>"),
	bindObj.pstvdFname = raw_input().strip("\n")

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
		for ci in range(nucLen):
			binaryNucleotide[4*ci:4*(ci + 1)] = map(int, encoding[nucleotide[ci]])
	else:
		for ci in range(nucLen):
			binaryNucleotide[4*((nucLen - ci) - 1): 4 * (nucLen - ci)] = map(int, reverseEncoding[nucleotide[ci]])

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
	with open(bindObj.pstvdFname, 'r') as pstvdFid:
		bindObj.pstvdSeq = pstvdFid.readline().strip("\n")
	
	bindObj.pstvdLen = len(bindObj.pstvdSeq)

	# Compute the total number of nucleotides and the maximum length of a nucleotide sequence
	nNucs = 0
	maxSeqLen = 0
	bindObj.nSequencesByLengths = np.zeros(bindObj.pstvdLen, dtype = int)
	with open(bindObj.sRNAPoolFname, 'r') as sRNAFid:
		for (lNo, line) in enumerate(sRNAFid):
			if ((lNo % 4) == 1):
				nNucs = nNucs + 1
				seqLen = len(line.strip("\n").strip(" "))

				bindObj.nSequencesByLengths[seqLen] = bindObj.nSequencesByLengths[seqLen] + 1

				if seqLen > maxSeqLen:
					maxSeqLen = seqLen

	bindObj.totalNucleotides = nNucs
	bindObj.maxNucleotideLength = maxSeqLen

	# Construct the relevant arrays
	bindObj.pstvdMatrix = np.zeros((4 * bindObj.maxNucleotideLength, bindObj.pstvdLen), dtype = int)
	
	if (((QUIET == 0)) or (DEBUG > 0)):
		print("\033[94mComputing the PSTVd encoding matrix.\033[0m")
	
	for si in range(bindObj.pstvdLen):
		
		if DEBUG > 1:
			print("\033[95mSlicing from %d to %d.\033[0m" % (si, si + bindObj.maxNucleotideLength))

		bindObj.pstvdMatrix[:, si] = NucleotideEncoder(StringSlice(bindObj.pstvdSeq, si, bindObj.maxNucleotideLength), bindObj.maxNucleotideLength, 0)
	
	if DEBUG > 1:
		print(bindObj.pstvdMatrix)
	
	bindObj.rawForwMatchFname = ("raw_forward_matchings_%s_%s_tol%d.txt" % (bindObj.sRNAPoolFname[:-4], bindObj.pstvdFname[:-4], bindObj.tolerance))
	with open(bindObj.rawForwMatchFname, 'w') as rawFid:
		rawFid.write("%d %d\n" % (bindObj.totalNucleotides, bindObj.pstvdLen))

	bindObj.rawRevMatchFname = ("raw_reverse_matchings_%s_%s_tol%d.txt" % (bindObj.sRNAPoolFname[:-4], bindObj.pstvdFname[:-4], bindObj.tolerance))
	with open(bindObj.rawForwMatchFname, 'w') as rawFid:
		rawFid.write("%d %d\n" % (bindObj.totalNucleotides, bindObj.pstvdLen))

	bindObj.forwardMatchFname = ("forward_matchings_%s_%s_tol%d.txt" % (bindObj.sRNAPoolFname[:-4], bindObj.pstvdFname[:-4], bindObj.tolerance))
	bindObj.forwardMatches = np.zeros((np.count_nonzero(bindObj.nSequencesByLengths), bindObj.pstvdLen + 1), dtype = int)
	bindObj.reverseMatchFname = ("reverse_matchings_%s_%s_tol%d.txt" % (bindObj.sRNAPoolFname[:-4], bindObj.pstvdFname[:-4], bindObj.tolerance))
	bindObj.reverseMatches = np.zeros((np.count_nonzero(bindObj.nSequencesByLengths), bindObj.pstvdLen + 1), dtype = int)

	bindObj.forwardMatches[:,0] = np.nonzero(bindObj.nSequencesByLengths)[0]
	bindObj.reverseMatches[:,0] = np.nonzero(bindObj.nSequencesByLengths)[0]

	bindObj.forwardMatchCounts = np.zeros(bindObj.pstvdLen, dtype = int)
	bindObj.reverseMatchCounts = np.zeros(bindObj.pstvdLen, dtype = int)
	return None


def StringSlice(seqStr, start, width):
	# Return the substring of a string that starts at a given index and has a given number of characters. The output can be a substring when the string is cyclic.
	substring = []
	for si in range(start, start + width):
		substring.append(seqStr[(si % len(seqStr))])
	return substring


def SetBreakPoints(bindObj):
	# Set points which will separate the nucleotides in the sRNA pool into chunks.
	# Each chunk will have at most MAXCONCURRENT nucleotides can be read by a core of the machine.
	nSeqsPartialPool = 0
	nSeqsRead = 0
	breakpoints = [[0,0,0]]
	with open(bindObj.sRNAPoolFname, 'r') as sRNAFid:
		for (lNo, line) in enumerate(sRNAFid):
			sRNASeq = line.strip("\n")
			# every 4th line contains a nucleotide sequence. Others contain data irrelevant to this algorithm
			if (lNo % 4 == 1):
				if ('N' in sRNASeq):
					continue
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

	return breakpoints


def PartialBinding(localMatchQueue, bindObj, nSkip, startSeq, numSeqs, reverse):
	# Partially load the sRNA pool with either the remaining or the next MAXCONCURRENT nucleotides, from the main pool, whichever is lesser.
	# The number of nucleotide sequences to be read can be less than the above estimate -- because valid nucleotide sequences should not have 'N' in them.
	sRNAPoolMatrix = np.zeros((numSeqs, 4 * bindObj.maxNucleotideLength), dtype = int)
	seqLenOffsets = np.zeros((numSeqs, bindObj.pstvdLen), dtype = int)
	sRNASeqLengths = np.zeros(numSeqs, dtype = int)
	
	directions = ["forward", "reverse"]

	if (((QUIET == 0)) or (DEBUG > 0)):
		print("\033[2mLoading the next %d sequences, starting at index %d, in the sRNA pool for %s pattern matching...\033[0m" % (numSeqs, startSeq, directions[reverse]))

	with open(bindObj.sRNAPoolFname, 'r') as sRNAFid:
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
					sRNASeq = line.strip("\n")
					if ('N' in sRNASeq):
						if DEBUG > 1:
							print("\033[94mSkipping sequence %s.\033[0m" % sRNASeq)
					else:
						if (nRead < numSeqs):
							if DEBUG > 1:
								print("\033[95mReading sequence %s.\033[0m" % sRNASeq)

							sRNAPoolMatrix[nRead, :] = NucleotideEncoder(StringSlice(sRNASeq, 0, len(sRNASeq)), bindObj.maxNucleotideLength, reverse)
							sRNASeqLengths[nRead] = len(sRNASeq)
							seqLenOffsets[nRead, :] = (len(sRNASeq) - bindObj.tolerance) * np.ones(bindObj.pstvdLen, dtype = int)
							nRead = nRead + 1

	if DEBUG > 1:
		print("sRNA subpool encoding matrix in the direction %d." % reverse)
		print(sRNAPoolMatrix.shape)
		print("sRNA sequence length matrix")
		print(sRNASeqLengths.shape)

	
	# Compute the frequency of matchings in the sRNA pool with every continuous PSTVd subsequence.
	# First, compute the Pattern matchings by matrix multiplication of the sRNA sequence encoding matrix with the PSTVd sequences encoding matrix.
	# If there was a perfect matching between a PSTVd subsequene and the sRNA nucleotide in the pool, the corresponding element generated in the product will be equal to the length of the sRNA nucleotide.
	# The amount by which the product element is lesser than the length of the sRNA nucleotide, is precisely the number of mismatches.
	# If the product doescribes matching upto the right tolerence, then reduce the product to 1, or else to 0.

	# matchings = matrix whose i,j entry is equal to the number of locations at which the ith nucleotide in the sRNA pool matched with the pstvd sequence of length L starting at index j.
	matchings = np.dot(sRNAPoolMatrix, bindObj.pstvdMatrix)
	# reducedToTolerance = matrix whose i,j entry is 1 if the ith nucleotide matched with the pstvd sequence (of length L) starting at index j.	
	reducedToTolerance = np.greater_equal(matchings, seqLenOffsets)

	(localNucIndices, pstvdIndices) = np.nonzero(reducedToTolerance)
	lengthCounts = np.bincount(sRNASeqLengths[localNucIndices], minlength = bindObj.pstvdLen)

	rawMatchingInfo = np.transpose(np.vstack((localNucIndices + startSeq * np.ones_like(localNucIndices), pstvdIndices)))
	if reverse == 0:
		WriteMatrixToFile(bindObj.rawForwMatchFname, rawMatchingInfo, 1, 0, 0)
	else:
		WriteMatrixToFile(bindObj.rawRevMatchFname, rawMatchingInfo, 1, 0, 0)

	localMatches = np.zeros((np.count_nonzero(bindObj.nSequencesByLengths), bindObj.pstvdLen), dtype = int)
	nzLen = 0
	
	for li in range(bindObj.pstvdLen):
		if (bindObj.nSequencesByLengths[li] > 0):
			localMatches[nzLen, :] = np.sum(reducedToTolerance * np.repeat(np.array(sRNASeqLengths == li, dtype = int)[:, np.newaxis], bindObj.pstvdLen, axis = 1), axis = 0)
			nzLen = nzLen + 1

	localMatchQueue.put([localMatches, lengthCounts])
	return None


def RecordBinding(bindObj):
	# Save the matching frequency array
	WriteMatrixToFile(bindObj.forwardMatchFname, bindObj.forwardMatches, 0, 0, 0)
	WriteMatrixToFile(bindObj.reverseMatchFname, bindObj.reverseMatches, 0, 0, 0)
	return None


def LogBinding(bindObj):
	# Store the details of the pattern matching run into a log file. The file is to be appened.
	logFile = "details.txt"
	with open(logFile, 'w') as logFid:
		logFid.write("\n***************************************\n")

		logFid.write("Log file for simulations done on, at, %s.\n" % bindObj.dateStamp)

		logFid.write("Length of the PSTVd sequence in %s: %d.\n" % (bindObj.pstvdFname, bindObj.pstvdLen))
		logFid.write("Number of nucleotides in the sRNA pool in %s: %d.\n" % (bindObj.sRNAPoolFname, bindObj.totalNucleotides))
		logFid.write("Number of mismatches allowed: %d.\n" % bindObj.tolerance)
		logFid.write("Length\tIn pool\tFoward\tReverse\n")
		for li in range(bindObj.pstvdLen):
			if (bindObj.nSequencesByLengths[li] > 0):
				logFid.write("%d\t%d\t%d\t%d\n" % (li, bindObj.nSequencesByLengths[li], bindObj.forwardMatchCounts[li], bindObj.reverseMatchCounts[li]))
		logFid.write("-------\n")
		logFid.write("Total:\t%d\t%d\t%d\n\n" % (bindObj.totalNucleotides, np.sum(bindObj.forwardMatchCounts), np.sum(bindObj.reverseMatchCounts)))

		# logFid.write("Forward matching data:\n")

	# WriteMatrixToFile(logFile, bindObj.forwardMatches, 1, 0, 0)
	
	# with open(logFile, 'a') as logFid:
		# logFid.write("Reverse matching data:\n")
	
	# WriteMatrixToFile(logFile, bindObj.reverseMatches, 1, 0, 0)

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
	rnas = Binder()
	# Reading the gene, pool files and tolerance.
	with open(sys.argv[1], "r") as fp:
		for (l, line) in enumerate(fp):
			if (l == int(sys.argv[2])):
				linecontents = map(lambda ln: ln.strip("\n").strip(" "), line.split())
				rnas.pstvdFname = linecontents[0]
				rnas.sRNAPoolFname = linecontents[1]
				rnas.tolerance = int(linecontents[2])
				nCores = int(linecontents[3])
	
	ComputeDerivedParameters(rnas)
	breakpoints = SetBreakPoints(rnas)
	nBreaks = len(breakpoints)

	# print("nBreaks = %d" % nBreaks)

	print("\033[93mMatching with %d sRNA nucleotides and a PSTVd sequence of length %d. Maximum length of a nucleotide sequence is %d.\033[0m" % (rnas.totalNucleotides, rnas.pstvdLen, rnas.maxNucleotideLength))
	startTime = time.time()
	
	for di in range(2):
		# change above 'range (2)' to have both forward and reverse matching
		if (di == 0):
			print("\033[94mForward matchings:")
		else:
			print("\033[94mReverse matchings:")

		completed = 0
		while (completed < nBreaks):
			iterationStartTime = time.time()
			launchNow = min(nBreaks - completed, nCores)
			localMatchQueue = Queue()
			processes = []

			if ((QUIET == 0)):
				print("\033[93mGoing to launch %d cores at a time.\033[0m" % launchNow)
				if (launchNow == 0):
					sys.stderr.write("nCores = 0\n")

			for pi in range(launchNow):
				linesToSkip = breakpoints[completed + pi][0]
				seqsReadSoFar = breakpoints[completed + pi][1]
				seqsToRead = breakpoints[completed + pi][2]
				processes.append(Process(target = PartialBinding, args = (localMatchQueue, rnas, linesToSkip, seqsReadSoFar, seqsToRead, di)))

			for pi in range(launchNow):
				processes[pi].start()

			# Gathering results
			for pi in range(launchNow):
				if (di == 0):
					localResults = localMatchQueue.get()
					rnas.forwardMatches[:,1:] = np.add(rnas.forwardMatches[:,1:], localResults[0])
					rnas.forwardMatchCounts = np.add(rnas.forwardMatchCounts, localResults[1])
				else:
					localResults = localMatchQueue.get()
					rnas.reverseMatches[:,1:] = np.add(rnas.reverseMatches[:,1:], localResults[0])
					rnas.reverseMatchCounts = np.add(rnas.reverseMatchCounts, localResults[1])
					pass

			for pi in range(launchNow):
				processes[pi].join()

			completed = completed + launchNow

			iterationEndTime = time.time()
			iterationTime = iterationEndTime - iterationStartTime
			# Predict the remaining time to finish all the matching: remaining = time taken for single iteration x number of iterations left
			# The time taken for a single iteration can be measured in the first iterationa and the number of iterations left is at most the total number of sequences in the pool divided by the number of sequences read in a chunk.
			remainingTime = iterationTime * (nBreaks - completed)

			if (QUIET == 0):
				print("\033[95m ......... %d%% done, about %d minutes %d seconds left\033[0m" % (100 * completed/float(nBreaks), (remainingTime/60), (int(remainingTime) % 60)))
			else:
				print("\r\033[95m ......... %d%% done, about %d minutes %d seconds left\033[0m" % (100 * completed/float(nBreaks), (remainingTime/60), (int(remainingTime) % 60))),
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
	for li in range(rnas.pstvdLen):
		if (rnas.nSequencesByLengths[li] > 0):
			print("\033[92m%d\t%d\t%d\t%d\033[0m" % (li, rnas.nSequencesByLengths[li], rnas.forwardMatchCounts[li], rnas.reverseMatchCounts[li]))
	print("\033[92m--------------\033[0m")
	print("\033[92mTotal:\t%d\t%d\t%d\033[0m" % (rnas.totalNucleotides, np.sum(rnas.forwardMatchCounts), np.sum(rnas.reverseMatchCounts)))
	
	# Save the matching information to a file
	RecordBinding(rnas)
	LogBinding(rnas)

	print("\033[93mSummary written to %s and %s.\033[0m" % (rnas.forwardMatchFname, rnas.reverseMatchFname))
