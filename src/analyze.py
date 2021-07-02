import os
from sys import stdout
import datetime as dt
from tqdm import tqdm
import numpy as np
try:
	import matplotlib
	matplotlib.use("Agg")
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.pyplot as plt
	matplotlib.rcParams['axes.linewidth'] = 1.5 # set the value globally
	matplotlib.rcParams['xtick.major.width'] = 1.5
	matplotlib.rcParams['ytick.major.width'] = 1.5
	matplotlib.rcParams['xtick.labelsize'] = 48
	matplotlib.rcParams['ytick.labelsize'] = 48
	matplotlib.rcParams["font.family"] = "arial"
except Exception:
	sys.stderr.write("\033[91mMatplotlib does not exist, cannot make plots.\n\033[0m")

DEBUG = 0

class Canvas():
	"""
	Contains all the information required for a plot

	name: name of the plot

	Filenames:
	forwMapDataFile: file containing data for forward mapping
	revMapDataFile: file containing data for reverse mapping
	forwMapPlotFile = plot for forward mapping
	revMapPlotFile = plot for reverse mapping

	scaledForwDataFile = scaled version of the reverse data file 
	scaledRevDataFile = scaled version of the reverse data file

	Plot preferences:
	title: title of the plot
	xLabel: label for x-axis
	yLabel: label for y-axis
	forwColor:
	revColor:

	integer variables:
	isSeparate = whether or not the plots for the reverse and forward mappings must be separated
	isReverse = whether or not the plots for the reverse mapping need to be produced.
	scaling = scaling of the y-axis to compare the plot with another.

	integer arrays
	forwardMatchData = data for the forward mapping
	reverseMatchData = data for the reverse mapping
	"""

	def __init__(self, gene, pool, tol):
		# with open(interfaceFile, 'r') as uiFid:
		# 1. file name of nucleotides in the pool
		# 2. file name of pstvd sequence
		# 3. tolerance
		# self.pool = ("./../data/input/%s" % pool)
		# self.gene = ("./../data/input/%s" % gene)
		# self.tol = int(tol)

		self.pool = pool
		self.gene = gene
		self.tol = int(tol)

		self.name = ("PSTVd sequence: %s and sRNA Pool used: %s." % (self.gene, self.pool))
		self.forwMatchPlotFile = ("./../plots/forward_matchings_%s_%s_tol%d.pdf" % (self.pool[:-4], self.gene[:-4], self.tol))
		self.revMatchPlotFile = ("./../plots/reverse_matchings_%s_%s_tol%d.pdf" % (self.pool[:-4], self.gene[:-4], self.tol))
		self.matchPlotFile = ("./../plots/matchings_%s_%s_tol%d.pdf" % (self.pool[:-4], self.gene[:-4], self.tol))

		self.forwMatchDataFile = ("./../data/output/forward_matchings_%s_%s_tol%d.txt" % (self.pool[:-4], self.gene[:-4], self.tol))
		self.revMatchDataFile = ("./../data/output/reverse_matchings_%s_%s_tol%d.txt" % (self.pool[:-4], self.gene[:-4], self.tol))

		# Scaled data files can be decided later, depending of the data repspect to which the scaling is to be done
		self.scaledForwDataFile = ("./../data/output/norm_forward_matchings_%s_%s_tol%d.txt" % (self.pool[:-4], self.gene[:-4], self.tol))
		self.scaledRevDataFile = ("./../data/output/norm_reverse_matchings_%s_%s_tol%d.txt" % (self.pool[:-4], self.gene[:-4], self.tol))
		self.scaledForwPlotFile = ("./../data/output/norm_forward_matchings_%s_%s_tol%d.pdf" % (self.pool[:-4], self.gene[:-4], self.tol))
		self.scaledRevPlotFile = ("./../data/output/norm_reverse_matchings_%s_%s_tol%d.pdf" % (self.pool[:-4], self.gene[:-4], self.tol))


		self.title = ("Mapping nucleotides in %s with PSTVd sequence in %s, allowing at most %d mismatches" % (self.pool, self.gene, self.tol))
		self.xLabel = ("Index on PSTVd sequence")
		self.yLabel = ("Number of bound nucleotides")
		self.forwMarker = "s"
		self.revMarker = "^"
		self.forwColor = "r"
		self.revColor = "b"
		self.linewidth = 0
		self.format = "pdf"

		self.scaling = 1
		self.isSeparate = 1
		self.isReverse = 0

		if os.path.isfile(self.forwMatchDataFile):
			self.forwardMatchData = ReadMatrixFromFile(self.forwMatchDataFile)
			self.gene_length = self.forwardMatchData.shape[1] - 1
		else:
			print("\033[95m!! Forward matching data file %s does not exist.\033[0m" % self.forwMatchDataFile)

		if os.path.isfile(self.scaledForwDataFile):
			self.nForw = ReadMatrixFromFile(self.scaledForwDataFile, dataType = 'f')
		else:
			self.nForw = np.zeros((self.forwardMatchData.shape[0], self.forwardMatchData.shape[1]), dtype = np.float)
			self.nForw = self.forwardMatchData
			# print("self.nForw\n%s" % np.array_str(self.nForw))

		if os.path.isfile(self.revMatchDataFile):
			self.reverseMatchData = ReadMatrixFromFile(self.revMatchDataFile)
		else:
			print("\033[95m!! Reverse matching data file %s does not exist.\033[0m" % self.revMatchDataFile)

		if os.path.isfile(self.scaledRevDataFile):
			self.nRev = ReadMatrixFromFile(self.scaledRevDataFile, dataType = 'f')
		else:
			self.nRev = np.zeros((self.reverseMatchData.shape[0], self.reverseMatchData.shape[1]), dtype = np.float)
			self.nRev = self.reverseMatchData
			# print("self.nForw\n%s" % np.array_str(self.nForw))
		# 
		# self.lengths = [[21], [22], [23], [24], [21, 22, 23, 24]]
		self.lengths = [[5], [6], [5, 6]] # Only for example
		self.xscale = np.array([0, self.gene_length, self.gene_length/5], dtype = int)
		self.yscale = np.array([[1, np.max(np.sum(self.forwardMatchData[:, 1:], axis = 0)), 10],
								[np.min(np.sum(self.reverseMatchData[:, 1:], axis = 0)), 10, -1]], dtype = int)

		# print "reverseMatchData"
		# print self.reverseMatchData
		# print "nRev"
		# print self.nRev[np.nonzero(self.nRev)]

		return None


def ScaleMatchings(plobj):
	# If there is a matching at position i, for length l nucleotide,
	# then we want to set: y[i + j] = max(y[i + j], y[i]), for all 0 < j < l.
	# print("Forward\n{}\nReverse\n{}".format(plobj.nForw, plobj.nRev))
	scaled_forward = np.zeros_like(plobj.nForw[:, 1:])
	scaled_reverse = np.zeros_like(plobj.nRev[:, 1:])
	for l in range(plobj.nForw.shape[0]):
		for i in range(plobj.gene_length):
			if (np.abs(plobj.nForw[l, 1 + i]) > 0):
				for j in range(int(plobj.nForw[l, 0])):
					scaled_forward[l, (i + j) % plobj.gene_length] = max(plobj.nForw[l, 1 + i], scaled_forward[l, (i + j) % plobj.gene_length])
					# scaled_forward[l, (i + j) % plobj.gene_length] += plobj.nForw[l, 1 + i]
			if (np.abs(plobj.nRev[l, 1 + i]) > 0):
				for j in range(int(plobj.nRev[l, 0])):
					scaled_reverse[l, (i + j) % plobj.gene_length] = max(plobj.nRev[l, 1 + i], scaled_reverse[l, (i + j) % plobj.gene_length])
					# scaled_reverse[l, (i + j) % plobj.gene_length] += plobj.nRev[l, 1 + i]
		# if (plobj.nForw[l, 0] > 0):
		# 	print("length: {}\nscaled_forward\n{}".format(plobj.nForw[l, 0], scaled_forward[l, 1:]))
		# if (plobj.nRev[l, 0] > 0):
		# 	print("length: {}\nscaled_reverse\n{}".format(plobj.nForw[l, 0], scaled_reverse[l, 1:]))
	scaled_forward = np.concatenate((plobj.nForw[:, 0, np.newaxis], scaled_forward), axis=1)
	scaled_reverse = np.concatenate((plobj.nRev[:, 0, np.newaxis], scaled_reverse), axis=1)
	return (scaled_forward, scaled_reverse)


def ReadMatrixFromFile(fname, dataType = 'i', atol = 10E-7):
	# Read a matrix from file
	# print("Going to read %s." % fname)
	with open(fname, 'r') as matFid:
		(nRows, nCols) = list(map(int, matFid.readline().strip("\n").split(" ")))
		mat = np.zeros((nRows, nCols), dtype = float)
		for ri in range(nRows):
			line = matFid.readline().strip("\n").strip(" ")
			lineContents = line.split(" ")
			if (dataType == 'i'):
				mat[ri,:] = list(map(int, lineContents))
			else:
				mat[ri,:] = list(map(float, lineContents))

	# Replace small values by 0
	mat[mat < atol] = 0
	return mat


def WriteMatrixToFile(fname, mat, is_append = 0, sparse = 0, binary = 0, dataType = 'i'):
	# Write a matrix to file
	(nRows, nCols) = mat.shape

	if is_append == 0:
		if os.path.isfile(fname):
			os.remove(fname)

	with open(fname, 'a') as matFid:
		if is_append == 0:
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
						if (dataType == 'i'):
							matFid.write("%d " % (mat[ri, ci]))
						else:
							matFid.write("%f " % (mat[ri, ci]))

			matFid.write("\n")

	return None


def Load(plobj):
	# Load the data required for plotting from a file
	plobj.forwardMatchData = ReadMatrixFromFile(plobj.forwMatchDataFile)
	plobj.reverseMatchData = ReadMatrixFromFile(plobj.revMatchDataFile)
	# print("Loading forward\n{}\nLoading reverse\n{}".format(plobj.forwardMatchData, plobj.reverseMatchData))
	return None


# def Save(plobj):
# 	# Save the (modified) data to a file
# 	WriteMatrixToFile(plobj.forwMatchDataFile, plobj.forwardMatchData)
# 	WriteMatrixToFile(plobj.revMatchDataFile, plobj.reverseMatchData)
# 	return None


def LeastGreatestMultiple(number, factor):
	# Given a number, n, compute the lowest multiple, c, of a factor, x, such that: n <= c.x.
	# The multiple c = ceil(n/x)
	if (np.abs(factor) >= np.abs(number)):
		return number
	if number < 0:
		multiple = np.floor(number/factor)
	else:
		multiple = np.ceil(number/factor)
	return (multiple * factor)

def IntArrayToString(arr):
	# Convert an integer array to string.
	return ",".join(list(map(lambda x: "%d" % x, arr)))


def Plot(plobj, isScaled = 0, quiet = 0):
	# Produce plots with the selected data files and the plot parameters
	# print("Forward\n{}\nReverse\n{}".format(plobj.nForw, plobj.nRev))
	alignments = ["forward", "reverse"]
	# fill_shade: RG: 0.3, I: 0.65, example: 0.65
	fill_shade = "0.65"
	xreps = 1
	shift = 0
	plotfname = ("./../plots/%s_%s_tol_%d.pdf" % (plobj.gene[:plobj.gene.index(".")], plobj.pool[:plobj.pool.index(".")], plobj.tol))
	xaxis = np.arange(plobj.gene_length + shift)
	xtick_step = 50
	# (scaled_forward, scaled_reverse) = (plobj.nForw, plobj.nRev)
	(scaled_forward, scaled_reverse) = ScaleMatchings(plobj)
	with PdfPages(plotfname) as pdf:
		for l in range(len(plobj.lengths)):
			for d in range(2):
				if (d == 0):
					nuclens = np.where(np.in1d(scaled_forward[:, 0], plobj.lengths[l]))[0]
					yaxis = np.tile(np.sum(scaled_forward[nuclens, 1:], axis = 0), xreps)[:len(xaxis)]
				else:
					nuclens = np.where(np.in1d(scaled_reverse[:, 0], plobj.lengths[l]))[0]
					yaxis = (-1) * np.tile(np.sum(scaled_reverse[nuclens, 1:], axis = 0), xreps)[:len(xaxis)]
				
				# Non-zero indices on the gene
				nonzero_gene_indicator = np.where(yaxis != 0, 1, 0)

				fig = plt.figure(figsize = (32, 24))
				# plt.plot(xaxis, yaxis, color = "0.5")
				currax = plt.gca()
				plt.fill_between(xaxis, np.zeros(yaxis.shape[0], dtype = np.float), yaxis, where=nonzero_gene_indicator, color=fill_shade)
				if (d == 0):
					# ######
					# Only for example
					currax.set_yticks(np.arange(0, 6, 1, dtype = np.int))
					currax.set_yticklabels(np.arange(0, 6, 1, dtype = np.int), fontsize = 54)
					currax.set_ylim([0, 5])
					# ######
					# y_top = max(6000, LeastGreatestMultiple(1.1 * np.max(yaxis), 500))
					# currax.set_yticks(np.arange(0, y_top, 1000))
					# currax.set_yticklabels(np.arange(0, y_top, 1000, dtype = np.int), fontsize = 54)
					# currax.set_ylim([0, 5000])
					# currax.set_ylim([0, LeastGreatestMultiple(1.05 * np.max(yaxis), 500)])
				else:
					currax.xaxis.tick_top()
					######
					# Only for example
					currax.set_yticks(np.arange(0, -6, -1, dtype = np.int))
					currax.set_yticklabels(np.arange(0, -6, -1, dtype = np.int), fontsize = 54)
					currax.set_ylim([-5, 0])
					######
					# y_bottom = min(-600, LeastGreatestMultiple(1.1 * np.max(yaxis), 100))
					# currax.set_yticks(np.arange(y_bottom, 0, 100, dtype = np.int))
					# currax.set_yticklabels(np.arange(y_bottom, 0, 100), fontsize = 54)
					# currax.set_ylim([-500, 0])
					# currax.set_ylim([LeastGreatestMultiple(1.05 * np.min(yaxis), 500), 0])
				
				# xticks and limits
				# xticks = list(xaxis[::50])
				# xticks.append(xaxis[-1])
				# currax.set_xticks(xticks)
				# currax.set_xticklabels([1 + xticks[0]] + xticks[1:-1] + [1 + xticks[-1]], rotation=45, fontsize = 54)
				# currax.set_xlim([0, None])
				#######
				# # For example
				currax.set_xticks(xaxis)
				currax.set_xticklabels(1 + np.array(xaxis, dtype = np.int), rotation=45, fontsize = 54)
				currax.set_xlim([0, len(xaxis) - 1])
				#######
				
				plt.title("%s matching for lengths %s" % (alignments[d], IntArrayToString(plobj.lengths[l])), fontsize = 48, y = 1.1)
				# print("xaxis\n{}\nnonzero\n{}".format(xaxis, np.nonzero(yaxis)))
				
				currax.tick_params(axis='y', which='both', pad = 20, direction = 'inout', length = 10, width = 3)
				# rightax.tick_params(axis='y', which='both', pad = 20, left = 'off', right = 'off', labeltop = 'off', labelbottom = 'off', labelright = 'off', labelleft = 'off', direction = 'in', length = 10, width = 3)
				# rightax.tick_params(axis='y', which='both', pad = 20, direction = 'in', length = 10, width = 3)
				if (d == 0):
					currax.tick_params(axis='x', which='both', pad = 20, direction = 'inout', length = 10, width = 3)
					# currax.tick_params(axis='x', which='both', pad = 20, direction = 'in', length = 10, width = 3)
				else:
					# currax.tick_params(axis='x', which='both', pad = 20, direction = 'inout', length = 10, width = 3)
					currax.tick_params(axis='x', which='both', pad = 20, direction = 'in', length = 10, width = 3)
				# currax.set_xticks(currax.get_xticks())
				# rightax.set_yticks(currax.get_yticks())
				
				print("Saving plot for length %s." % (ArrayToString(plobj.lengths[l])))
				pdf.savefig(fig)  # saves the current figure into a pdf page
				plt.close()
				# topax.cla()
				# rightax.cla()
		# Set the PDF attributes
		pdfInfo = pdf.infodict()
		pdfInfo['Title'] = ("Comparing different sampling methods")
		pdfInfo['Author'] = "Pavithran Iyer"
		pdfInfo['ModDate'] = dt.datetime.today()
	print("\033[92mAll plots saved to %s.\033[0m" % plotfname)
	return None


def ArrayToString(arr, data_type = 'i'):
	# Convert an array to string.
	if (data_type == 'i'):
		arr_str = ",".join(list(map(lambda x: "%d" % x, arr)))
	return arr_str


def Normalize(dataset):
	# Multiply the number of reads at every position by a constant that depends on the nucleatide length.
	# This is to ensure that all reads are scaled for a pool size that is 1,000,000.
	# Scaling array provides the scaling constants by pool gene tol length.
	million = 1E6
	pool_size = {
					"sRNA_intermediate.txt": {21: 700131, 22: 839033, 23: 622373, 24: 2155006},
					"GSM1717894_PSTVd_RG1.txt": {21: 160262, 22: 145487, 23: 124625, 24: 300125}
				}

	# Scale the forward matchings
	scaled_forward = np.zeros_like(dataset.forwardMatchData)
	scaled_forward[:, 0] = dataset.forwardMatchData[:, 0]
	for s in range(len(dataset.lengths)):
		nuclens = np.where(np.in1d(dataset.nForw[:, 0], dataset.lengths[s]))[0]
		total_seqs = np.sum([pool_size[dataset.pool][l] for l in dataset.lengths[s]])
		scaling = million/total_seqs
		scaled_forward[nuclens, 1:] = dataset.forwardMatchData[nuclens, 1:] * scaling
	
	WriteMatrixToFile(dataset.scaledForwDataFile, scaled_forward, is_append = 0, dataType='f')
	
	# Scale the forward matchings
	scaled_reverse = np.zeros_like(dataset.reverseMatchData)
	scaled_reverse[:, 0] = dataset.reverseMatchData[:, 0]
	for s in range(len(dataset.lengths)):
		nuclens = np.where(np.in1d(dataset.nForw[:, 0], dataset.lengths[s]))[0]
		total_seqs = np.sum([pool_size[dataset.pool][l] for l in dataset.lengths[s]])
		scaling = million/total_seqs
		scaled_reverse[nuclens, 1:] = dataset.reverseMatchData[nuclens, 1:] * scaling

	WriteMatrixToFile(dataset.scaledRevDataFile, scaled_reverse, is_append = 0, dataType='f')
	return None


def ListForwardMatching(dset):
	# List all the forward matching sequences
	# If the forward matching array element, F[i][j] = x, then we need the gene-substring gene[i:(i + lengths[i])] and x.
	with open("./../data/input/%s" % (dset.gene), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	forward_matches_log = "./../data/output/explicit_forward_%s_%s_%d.txt" % (dset.gene, dset.pool, dset.tol)
	with open(forward_matches_log, "w") as fl:
		fl.write("Forward matchings\n\n")
		fl.write("Gene: %s\n" % (dset.gene))
		fl.write("Pool: %s\n" % (dset.pool))
		fl.write("Mismatches: %d\n" % (dset.tol))
		fl.write("*************************\n\n")
		for l in range(dset.forwardMatchData.shape[0]):
			nuc_len = int(dset.forwardMatchData[l, 0])
			gene_indices, = np.nonzero(dset.forwardMatchData[l, 1:])
			match_freq = dset.forwardMatchData[l, 1 + gene_indices].astype(np.int)
			if (gene_indices.shape[0] > 0):
				fl.write("Length: %d\n" % (nuc_len))
				fl.write("{:^12} | {:^8}\n".format("Sequence", "Frequency"))
				fl.write("-------------------------\n")
				for s in range(gene_indices.shape[0]):
					gene_subseq = [gene_seq[(gene_indices[s] + g) % len(gene_seq)] for g in range(nuc_len)]
					fl.write("{:^12} | {:^8}\n".format("".join(gene_subseq), match_freq[s]))
				fl.write("*************************\n\n")
	return None


def ReverseComplement(seq):
	# Compute the reverse complement encoding.
	reverse_encoding = {"A":"T", "T":"A", "G":"C", "C":"G"}
	revcomp = [seq[s] for s in range(len(seq))]
	for s in range(len(seq)):
		revcomp[s] = reverse_encoding[seq[s]]
	return "".join(revcomp[::-1])


def ListReverseMatching(dset):
	# List all the reverse matching sequences
	# If the reverse matching array element, F[i][j] = x, then we need the gene-substring gene[i:(i + lengths[i])] and x.
	with open("./../data/input/%s" % (dset.gene), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	reverse_matches_log = "./../data/output/explicit_reverse_%s_%s_%d.txt" % (dset.gene, dset.pool, dset.tol)
	with open(reverse_matches_log, "w") as fl:
		fl.write("Reverse matchings\n\n")
		fl.write("Gene: %s\n" % (dset.gene))
		fl.write("Pool: %s\n" % (dset.pool))
		fl.write("Mismatches: %d\n" % (dset.tol))
		fl.write("*************************\n\n")
		for l in range(dset.reverseMatchData.shape[0]):
			nuc_len = int(dset.reverseMatchData[l, 0])
			gene_indices, = np.nonzero(dset.reverseMatchData[l, 1:])
			match_freq = dset.reverseMatchData[l, 1 + gene_indices].astype(np.int)
			if (gene_indices.shape[0] > 0):
				fl.write("Length: %d\n" % (nuc_len))
				fl.write("{:^12} | {:^8}\n".format("Sequence", "Frequency"))
				fl.write("-------------------------\n")
				for s in range(gene_indices.shape[0]):
					gene_subseq = ReverseComplement([gene_seq[(gene_indices[s] + g) % len(gene_seq)] for g in range(nuc_len)])
					fl.write("{:^12} | {:^8}\n".format("".join(gene_subseq), match_freq[s]))
				fl.write("*************************\n\n")
	return None


def ListMatching(dset):
	# List the forward and reverse matching output.
	ListForwardMatching(dset)
	ListReverseMatching(dset)
	return None


def VerifyForwardMatching(dset):
	# Verify the matching output from the dataset.
	# If the forward matching array element, F[i][j] = x, then we need to check if indeed the gene-substring gene[i:(i + lengths[i])] occurs x times in the pool.
	with open("./../data/input/%s" % (dset.gene), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	for l in range(len(dset.lengths)):
		is_valid = 1
		nuc_len = dset.lengths[l][0]
		for i in tqdm(range(1, dset.forwardMatchData.shape[1]), desc="Forward matchings for length %d" % (nuc_len)):
			nuc_len_idx = np.where(np.in1d(dset.forwardMatchData[:, 0], dset.lengths[l]))[0]
			gene_subseq = gene_seq[(i - 1)  : (i - 1 + nuc_len) % len(gene_seq)]
			given_matches = dset.forwardMatchData[nuc_len_idx, i]
			found_matches = MatchesInFile("./../data/input/%s" % (dset.pool), gene_subseq)
			
			if given_matches == found_matches:
				is_match = 1
				# print("[_/] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			else:
				is_match = 0
				print("[X] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			
			is_valid *= is_match
		
		if (is_valid == 1):
			print("vbind gave correct forward matching output for pool sequences with {} nucleotides.".format(nuc_len))
		else:
			print("vbind gave incorrect forward matching output for pool sequences with {} nucleotides.".format(nuc_len))
	return is_valid


def VerifyReverseMatching(dset):
	# Verify the reverse matching output
	# If the reverse matching array element, R[i][j] = x, then we need to check if indeed the reverse complement of the gene-substring gene[i:(i + lengths[i])] occurs x times in the pool.
	with open("./../data/input/%s" % (dset.gene), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	for l in range(len(dset.lengths)):
		is_valid = 1
		nuc_len = dset.lengths[l][0]
		for i in tqdm(range(1, dset.reverseMatchData.shape[1]), desc="Reverse matchings for length %d" % (nuc_len)):
			nuc_len_idx = np.where(np.in1d(dset.reverseMatchData[:, 0], dset.lengths[l]))[0]
			gene_subseq = ReverseComplement(gene_seq[(i - 1)  : (i - 1 + nuc_len) % len(gene_seq)])
			given_matches = dset.reverseMatchData[nuc_len_idx, i]
			found_matches = MatchesInFile("./../data/input/%s" % (dset.pool), gene_subseq)
			
			if given_matches == found_matches:
				is_match = 1
				# print("[_/] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			else:
				is_match = 0
				print("[X] Matches for {} of length {} at position {} on the gene.\nReported: {} and Found: {}.".format(gene_subseq, nuc_len, l, given_matches, found_matches))
			
			is_valid *= is_match
		
		if (is_valid == 1):
			print("vbind gave correct reverse matching output for pool sequences with {} nucleotides.")
		else:
			print("vbind gave incorrect reverse matching output for pool sequences with {} nucleotides.")
	return is_valid


def ReverseComplement(nuc_seq):
	# Compute the reverse complement of a nucleotide sequence
	mapping = {"A":"T", "T":"A", "G":"C", "C":"G"}
	complement = []
	for char in nuc_seq[::-1]:
		complement.append(mapping[char])
	return "".join(complement)


def VerifyMatching(dset):
	# Verify the forward and reverse matching output.
	is_valid = VerifyForwardMatching(dset) * VerifyReverseMatching(dset)
	return is_valid


def MatchesInFile(fname, search_string):
	# Find the number of lines in a file which contain a given string.
	occurrences = 0
	# print("Searching for a matching for {} in the file {}".format(search_string, fname))
	with open(fname, "r") as f:
		for line in f:
			if (line.strip(" ").strip("\n") == search_string):
				occurrences += 1
	return occurrences


if __name__ == '__main__':
	completed = 0
	isMenu = 1
	canvanses = {}
	alive_time = 0
	# list of plot data for which normalized versions are available
	scaledAvailable = []
	while completed == 0:
		if (isMenu == 1):
			print("\033[92m**** MENU ****\033[0m")
			print("\033[93m0 -- Quit\033[0m")
			print("\033[93m1 -- Load new data for matching\033[0m")
			print("\033[93m2 -- Plot a previously loaded data\033[0m")
			print("\033[93m3 -- Relative scaling of matching data\033[0m")
			print("\033[93m4 -- Write matching summary to a file\033[0m")
			print("\033[93m5 -- Verify matching")
			print("\033[93m6 -- Menu")
			isMenu = 0
		else:
			print("\033[2m---enter 6 for the menu---\033[0m")

		print(">>"),
		user_choice = int(input().strip("\n").strip(" "))

		# ####### Temporary
		# if alive_time == 0:
		# 	user_choice = 1
		# elif alive_time == 1:
		# 	user_choice = 4
		# elif alive_time == 2:
		# 	user_choice = 2
		# else:
		# 	user_choice = 0
		# #######

		if (user_choice == 0):
			completed = 1

		elif (user_choice == 1):
			# I & RG: +: 5000, tick at 1000
			# 		-: -500, tick at 100
			# inputs = [("pstvdI.txt", "sRNA_intermediate.txt", 0)]
			# inputs = [("pstvdI.txt", "sRNA_intermediate.txt", 1)]
			# inputs = [("PSTVd_RG.txt", "GSM1717894_PSTVd_RG1.txt", 0)]
			# inputs = [("PSTVd_RG.txt", "GSM1717894_PSTVd_RG1.txt", 1)]
			# inputs = [("example_gene.txt", "example_pool.txt", 0)]
			inputs = [("example_gene.txt", "example_pool.txt", 1)]
			name = "default"
			for (gene, pool, tol) in inputs:
				newCan = Canvas(gene, pool, tol)
				Load(newCan)
				canvanses.update({name:newCan})

		elif (user_choice == 2):
			print("\033[93mOf the %d available matching data sets, which one would you like to plot?\033[0m" % len(canvanses.keys()))
			scaled = 0
			for (ni, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (ni, name)),
				if (name in scaledAvailable):
					print("\033[92m\t (*)\033[0m")
					scaled = 1
				else:
					print("")
			print("\033[92m* -- Normalized forms available.\033[0m")
			
			print(">>"),
			# userSelection = input().strip("\n").strip(" ")
			userSelection = 0 # Temporary, remove finally.

			plotDataChoice = [int(userSelection)]
			
			for pdi in plotDataChoice:
				Plot(canvanses[list(canvanses.keys())[pdi]], isScaled = 1, quiet = 1)

		elif (user_choice == 3):
			print("\033[93mOf the %d available matching data sets, which ones do you want to include?\033[0m" % len(canvanses.keys()))
			for (ni, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (ni, name))
			print(">>"),
			# userSelection = int(input().strip("\n").strip(" "))
			userSelection = 0 ### temporary
			
			dataset = canvanses[list(canvanses.keys())[userSelection]]
			if os.path.isfile(dataset.scaledForwDataFile):
				os.remove(dataset.scaledForwDataFile)
			if os.path.isfile(dataset.scaledRevDataFile):
				os.remove(dataset.scaledRevDataFile)
			Normalize(dataset)

		elif (user_choice == 4):
			ListMatching(canvanses["default"])

		elif (user_choice == 5):
			print("\033[93mOf the %d available matching data sets, for which one would you like to a report?\033[0m" % len(canvanses.keys()))
			for (ni, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (ni, name))
			
			print(">>"),
			reportDataChoice = int(input().strip("\n").strip(" "))
			Report(canvanses[list(canvanses.keys())[reportDataChoice]])

		elif (user_choice == 7):
			isMenu = 1

		elif (user_choice == 6):
			print("\033[93mOf the %d available matching data sets, for which one would you like to verify?\033[0m" % len(canvanses.keys()))
			for (n, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (n, name))
			
			# print(">>"),
			# dset_choice = int(input().strip("\n").strip(" "))
			# Temporary
			dset_choice = 0
			VerifyMatching(canvanses[list(canvanses.keys())[dset_choice]])

		else:
			print("\033[91mUnknown option!\033[0m")

		alive_time += 1

	print("\033[92mxxxxxxxx\033[0m")
