import os
import sys
import numpy as np
from binder import GroupByLength
from normalize import RecordNormalizedData, Normalize
from verify import VerifyMatching
from utils import ArrayToString, IntArrayToString, ReverseComplement, LeastGreatestMultiple, ReadMatrixFromFile
from plot import PlotSettings, Plot
from startup import DisplayLogoLicense, CheckDependencies

class RNAPostProcessor():
	"""
	Contains all the information required for post processing the output of the sRNA profiler.
	"""

	def __init__(self):
		# Basic details
		self.poolfname = None
		self.genefname = None
		self.tolerance = 0
		self.file_skip = 0
		self.is_circular = 0
		self.poolsize = 0
		self.max_nuc_length = 0
		self.nSequencesByLengths = None
		self.scaled_forward = None
		self.scaled_reverse = None
		self.lengths = None
		return None

	def load_output(self):
		# Set the filenames for the matching data and load data
		# Matching data.
		self.forwMatchDataFile = ("./../data/output/forward_matchings_%s_%s_tol%d.txt" % (self.poolfname[:-4], self.genefname[:-4], self.tolerance))
		self.nForw = ReadMatrixFromFile(self.forwMatchDataFile, dataType = 'f')
		self.gene_length = self.nForw.shape[1] - 1

		self.revMatchDataFile = ("./../data/output/reverse_matchings_%s_%s_tol%d.txt" % (self.poolfname[:-4], self.genefname[:-4], self.tolerance))
		self.nRev = ReadMatrixFromFile(self.revMatchDataFile, dataType = 'f')
		
		# Normalized matching data.
		self.normalized_forward = None
		self.norm_ForwDataFile = ("./../data/output/norm_forward_matchings_%s_%s_tol%d.txt" % (self.poolfname[:-4], self.genefname[:-4], self.tolerance))
		
		self.normalized_reverse = None
		self.norm_RevDataFile = ("./../data/output/norm_reverse_matchings_%s_%s_tol%d.txt" % (self.poolfname[:-4], self.genefname[:-4], self.tolerance))
		return None


def Usage():
	# Print the usage.
	print("\033[2m./analyze.py <input file>\033[0m")
	print("\033[2mwhere the input file should be in data/input/ and the same format as for vbind.sh.\033[0m")
	return None


def IdentifyInstance(rnap, input_file, line_number):
	# Extract the parameters from an input file that specify an instance for post-processing.
	with open(input_file, "r") as fp:
		lno = 1
		for line in fp:
			if (line[0] != "#"):
				if (lno == line_number):
					line = list(map(lambda ln: ln.strip("\n").strip(" "), line.split(" ")))
					rnap.genefname = line[0]
					rnap.poolfname = line[1]
					rnap.file_skip = int(line[2])
					rnap.tolerance = int(line[3])
					rnap.is_circular = int(line[4])
					ncores = int(line[5])
					break
				lno = lno + 1
	rnap.load_output()
	return rnap


def ScaleMatchings(plobj, forward, reverse):
	# If there is a matching at position i, for length l nucleotide,
	# then we want to set: y[i + j] = max(y[i + j], y[i]), for all 0 < j < l.
	# print("plobj.gene_length = {}\nforward: shape = {}\n{}\nreverse: shape = {}\n{}".format(plobj.gene_length, forward.shape, forward, reverse.shape, reverse))
	plobj.scaled_forward = np.zeros_like(forward)
	plobj.scaled_reverse = np.zeros_like(reverse)
	for l in range(len(plobj.lengths)):
		for i in range(plobj.gene_length):
			if (np.abs(forward[l, i]) > 0):
				for j in range(sum(plobj.lengths[l])):
					plobj.scaled_forward[l, (i + j) % plobj.gene_length] = max(forward[l, i], plobj.scaled_forward[l, (i + j) % plobj.gene_length])
			if (np.abs(reverse[l, i]) > 0):
				for j in range(sum(plobj.lengths[l])):
					plobj.scaled_reverse[l, (i + j) % plobj.gene_length] = min(reverse[l, i], plobj.scaled_reverse[l, (i + j) % plobj.gene_length])
	return None


def Load(plobj):
	# Load the data required for plotting from a file
	plobj.forwardMatchData = ReadMatrixFromFile(plobj.forwMatchDataFile)
	plobj.reverseMatchData = ReadMatrixFromFile(plobj.revMatchDataFile)
	return None


def GatherMatchingData(rnap, forward, reverse):
	# Compute the number of matchings per each length set.
	rnap.gathered_forward = np.zeros((len(rnap.lengths), forward.shape[1] - 1), dtype = np.int)
	rnap.gathered_reverse = np.zeros((len(rnap.lengths), reverse.shape[1] - 1), dtype = np.int)
	for s in range(len(rnap.lengths)):
		for d in range(2):
			if (d == 0):
				nuclens = np.where(np.in1d(forward[:, 0], rnap.lengths[s]))[0]
				rnap.gathered_forward[s, :] = np.sum(forward[nuclens, 1:], axis = 0)
			else:
				nuclens = np.where(np.in1d(reverse[:, 0], rnap.lengths[s]))[0]
				rnap.gathered_reverse[s, :] = (-1) * np.sum(reverse[nuclens, 1:], axis = 0)
	return None


def SummarizeForwardMatching(dset):
	# List all the forward matching sequences
	# If the forward matching array element, F[i][j] = x, then we need the gene-substring gene[i:(i + lengths[i])] and x.
	with open("./../data/input/%s" % (dset.genefname), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	forward_matches_log = "./../data/output/explicit_forward_%s_%s_%d.txt" % (dset.genefname, dset.poolfname, dset.tolerance)
	topology = ["linear", "circular"][dset.is_circular]
	with open(forward_matches_log, "w") as fl:
		fl.write("Forward matchings\n\n")
		fl.write("Gene: %s\n" % (dset.genefname))
		fl.write("Pool: %s\n" % (dset.poolfname))
		fl.write("Topology: %s\n" % (topology))
		fl.write("Mismatches: %d\n" % (dset.tolerance))
		fl.write("*************************\n\n")
		for l in range(dset.nForw.shape[0]):
			nuc_len = int(dset.nForw[l, 0])
			gene_indices, = np.nonzero(dset.nForw[l, 1:])
			match_freq = dset.nForw[l, 1 + gene_indices].astype(np.int)
			if (gene_indices.shape[0] > 0):
				fl.write("Length: %d\n" % (nuc_len))
				fl.write("{:^12} | {:^8}\n".format("Sequence", "Frequency"))
				fl.write("-------------------------\n")
				for s in range(gene_indices.shape[0]):
					gene_subseq = [gene_seq[(gene_indices[s] + g) % len(gene_seq)] for g in range(nuc_len)]
					fl.write("{:^12} | {:^8}\n".format("".join(gene_subseq), match_freq[s]))
				fl.write("*************************\n\n")
	return None


def SummarizeReverseMatching(dset):
	# List all the reverse matching sequences
	# If the reverse matching array element, F[i][j] = x, then we need the gene-substring gene[i:(i + lengths[i])] and x.
	with open("./../data/input/%s" % (dset.genefname), "r") as gf:
		gene_seq = gf.readline().strip(" ").strip("\n")

	reverse_matches_log = "./../data/output/explicit_reverse_%s_%s_%d.txt" % (dset.genefname, dset.poolfname, dset.tolerance)
	topology = ["linear", "circular"][dset.is_circular]
	with open(reverse_matches_log, "w") as fl:
		fl.write("Reverse matchings\n\n")
		fl.write("Gene: %s\n" % (dset.genefname))
		fl.write("Pool: %s\n" % (dset.poolfname))
		fl.write("Topology: %s\n" % (topology))
		fl.write("Mismatches: %d\n" % (dset.tolerance))
		fl.write("*************************\n\n")
		for l in range(dset.nRev.shape[0]):
			nuc_len = int(dset.nRev[l, 0])
			gene_indices, = np.nonzero(dset.nRev[l, 1:])
			match_freq = dset.nRev[l, 1 + gene_indices].astype(np.int)
			if (gene_indices.shape[0] > 0):
				fl.write("Length: %d\n" % (nuc_len))
				fl.write("{:^12} | {:^8}\n".format("Sequence", "Frequency"))
				fl.write("-------------------------\n")
				for s in range(gene_indices.shape[0]):
					gene_subseq = ReverseComplement([gene_seq[(gene_indices[s] + g) % len(gene_seq)] for g in range(nuc_len)])
					fl.write("{:^12} | {:^8}\n".format("".join(gene_subseq), match_freq[s]))
				fl.write("*************************\n\n")
	return None


def SummarizeMatching(dset):
	# List the forward and reverse matching output.
	SummarizeForwardMatching(dset)
	SummarizeReverseMatching(dset)
	return None


def ParseNucLengths(lengths_encoding):
	# Parse the string input specifying the nucleotide lengths.
	# The nucleotide lengths is a list of lists.
	# Each list in the string is separated by a semicolon ";" and each element of a list is separated by a comma ",".
	lengths_string = list(map(lambda ln: ln.strip("\n").strip(" ").split(","), lengths_encoding.strip("\n").strip(" ").split(";")))
	lengths = [list(map(int, ln)) for ln in lengths_string]
	return lengths


if __name__ == '__main__':
	
	# Display the logo and license information
	DisplayLogoLicense()
	# Check if all the required packages exist
	CheckDependencies()

	rnap = RNAPostProcessor()
	# Read the parameters to identify the instance for post-processing
	if (len(sys.argv) < 2):
		Usage()
		exit(0)
	else:
		input_file = ("./../data/input/%s" % sys.argv[1].strip("\n"))

	completed = 0
	user_choice = 6
	while (completed == 0):
		if (user_choice == 1):
			# inputs = [("example_gene.txt", "example_pool.txt", 1)]
			instance = int(input(">>Problem instance from the input file %s: " % (os.path.basename(input_file))).strip("\n").strip(" "))
			IdentifyInstance(rnap, input_file, instance)
			rnap.pool_lengths = GroupByLength(rnap)

		elif (user_choice == 2):
			rnap.lengths = ParseNucLengths(input(">>Lengths: ").strip("\n").strip(" "))
			# print("lengths: {}".format(rnap.lengths))
			GatherMatchingData(rnap, rnap.nForw, rnap.nRev)
			is_scaled = int(input(">>Plot normalized data? [1]Yes, [0]No: ").strip("\n").strip(" "))
			if (is_scaled == 1):
				Normalize(rnap, rnap.gathered_forward, rnap.gathered_reverse)
				# GatherMatchingData(rnap, rnap.gathered_forward, rnap.gathered_reverse)
				(rnap.gathered_forward, rnap.gathered_reverse) = (rnap.normalized_forward, rnap.normalized_reverse)
			ScaleMatchings(rnap, rnap.gathered_forward, rnap.gathered_reverse)
			# print("Reverse\n{}".format(rnap.gathered_reverse))
			# Load plot settings
			plot_settings = PlotSettings()
			settings_fname = input(">>Settings file name (leave blank for default): ").strip("\n").strip(" ")
			if (len(settings_fname) == 0):
				settings_fname = "./../data/input/default_plot_settings.txt"
			plot_settings.load(settings_fname)
			Plot(rnap, plot_settings)

		elif (user_choice == 3):
			rnap.lengths = ParseNucLengths(input(">>Lengths: ").strip("\n").strip(" "))
			RecordNormalizedData(rnap)

		elif (user_choice == 4):
			topology = ["linear", "circular"][rnap.is_circular]
			print("Gene: %s" % (rnap.genefname))
			print("Pool: %s" % (rnap.poolfname))
			print("Topology: %s" % (topology))
			print("Mismatches: %d" % (rnap.tolerance))
			SummarizeMatching(rnap)

		elif (user_choice == 5):
			rnap.lengths = ParseNucLengths(input(">>Lengths: ").strip("\n").strip(" "))
			VerifyMatching(rnap)

		elif (user_choice == 6):
			print("**** MENU ****")
			print("0 -- Quit")
			print("1 -- Load new data for matching.")
			print("2 -- Plot the latest dataset.")
			print("3 -- Normalize matching data.")
			print("4 -- Save the matching summary to a file.")
			print("5 -- Verify matching results.")
			print("6 -- Show menu")
			print("**** MENU ****")

		else:
			pass

		print("\033[2m---Enter 6 to show the menu---\033[0m")

		user_input = input(">>Menu Option: ").strip("\n").strip(" ")
		if user_input.isnumeric():
			user_choice = int(user_input)
		else:
			user_choice = -1

		if (user_choice == 0):
			completed = 1
		

	print("\033[2mxxxxxxxx\033[0m")
