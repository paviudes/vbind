import os
from sys import stdout
import datetime as dt
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
		# self.sRNAPoolFname = ("./../data/input/%s" % pool)
		# self.pstvdFname = ("./../data/input/%s" % gene)
		# self.tolerance = int(tol)

		self.sRNAPoolFname = pool
		self.pstvdFname = gene
		self.tolerance = int(tol)

		self.name = ("PSTVd sequence: %s and sRNA Pool used: %s." % (self.pstvdFname, self.sRNAPoolFname))
		self.forwMatchPlotFile = ("./../plots/forward_matchings_%s_%s_tol%d.pdf" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))
		self.revMatchPlotFile = ("./../plots/reverse_matchings_%s_%s_tol%d.pdf" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))
		self.matchPlotFile = ("./../plots/matchings_%s_%s_tol%d.pdf" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))

		self.forwMatchDataFile = ("./../data/output/forward_matchings_%s_%s_tol%d.txt" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))
		self.revMatchDataFile = ("./../data/output/reverse_matchings_%s_%s_tol%d.txt" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))

		# Scaled data files can be decided later, depending of the data repspect to which the scaling is to be done
		self.scaledForwDataFile = ("./../data/output/norm_forward_matchings_%s_%s_tol%d.txt" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))
		self.scaledRevDataFile = ("./../data/output/norm_reverse_matchings_%s_%s_tol%d.txt" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))
		self.scaledForwPlotFile = ("./../data/output/norm_forward_matchings_%s_%s_tol%d.pdf" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))
		self.scaledRevPlotFile = ("./../data/output/norm_reverse_matchings_%s_%s_tol%d.pdf" % (self.sRNAPoolFname[:-4], self.pstvdFname[:-4], self.tolerance))


		self.title = ("Mapping nucleotides in %s with PSTVd sequence in %s, allowing at most %d mismatches" % (self.sRNAPoolFname, self.pstvdFname, self.tolerance))
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
			self.pstvdLength = self.forwardMatchData.shape[1] - 1
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

		self.lengths = np.array([[3]], dtype = int)
		self.xscale = np.array([0, self.pstvdLength, self.pstvdLength/5], dtype = int)
		self.yscale = np.array([[1, np.max(np.sum(self.forwardMatchData[:, 1:], axis = 0)), 10],
								[np.min(np.sum(self.reverseMatchData[:, 1:], axis = 0)), 10, -1]], dtype = int)

		# print "reverseMatchData"
		# print self.reverseMatchData
		# print "nRev"
		# print self.nRev[np.nonzero(self.nRev)]

		return None

def ReadMatrixFromFile(fname, dataType = 'i', atol = 10E-7):
	# Read a matrix from file
	# print("Going to read %s." % fname)
	with open(fname, 'r') as matFid:
		(nRows, nCols) = map(int, matFid.readline().strip("\n").split(" "))
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


def WriteMatrixToFile(fname, mat, appendMode, sparse = 0, binary = 0, dataType = 'i'):
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
	print("Loading forward\n{}\nLoading reverse\n{}".format(plobj.forwardMatchData, plobj.reverseMatchData))
	return None

def Save(plobj):
	# Save the (modified) data to a file
	WriteMatrixToFile(plobj.forwMatchDataFile, plobj.forwardMatchData)
	WriteMatrixToFile(plobj.revMatchDataFile, plobj.reverseMatchData)
	return None


def PromptAdvancedParameters(plobj):
	# Prompt the user for selecting advanced parameters in the plot. There are parameters that typically need not be modified.
	availableFormats = ["pdf", "jpg", "jpeg", "png"]

	print("\033[2m###############\033[0m")
	print("\033[93m0. Don't modify anything\033[0m")
	print("\033[93m1. Separate forward and reverse plots: \033[92m%d.\033[0m" % plobj.isSeparate)
	print("\033[93m2. Plot reverse mapping? \033[92m%d.\033[0m" % plobj.isReverse)
	print("\033[93m3. Forward area shade: \033[92m%s.\033[0m" % plobj.forwColor)
	print("\033[93m4. Forward marker \033[92m%s.\033[0m" % plobj.forwMarker)
	print("\033[93m5. Reverse area shade: \033[92m%s.\033[0m" % plobj.revColor)
	print("\033[93m6. Reverse marker \033[92m%s.\033[0m" % plobj.revMarker)
	print("\033[93m7. Line width = \033[92m%d.\033[0m" % plobj.linewidth)
	print("\033[93m8. Save plots in the format: \033[92m%s.\033[0m" % plobj.format)

	while completed == 0:
		print("\033[2mEnter "'"option number;new value"'"\033[0m")
		print(">>"),
		stdout.flush()
		choices = input().strip("\n").strip(" ")
		if len(choices) > 1:
			(optionStr, newValue) = choices.split(";")
			option = int(optionStr.strip(" "))
		else:
			exitOption = int(choices.strip(" "))
			print("\033[2mParameters remain in thier current state.\033[0m")
			option = 0

		if (option == 0):
			completed = 1

		elif (option == 1):
			plobj.isSeparate = int(newValue.strip(" "))
		elif (option == 2):
			plobj.isReverse = int(newValue.strip(" "))
		elif (option == 3):
			plobj.forwColor = newValue.strip(" ")
		elif (option == 4):
			plobj.forwMarker = newValue.strip(" ")
		elif (option == 5):
			plobj.revColor = newValue.strip(" ")
		elif (option == 6):
			plobj.revMarker = newValue.strip(" ")
		elif (option == 7):
			plobj.linewidth = int(newValue.strip(" "))
		elif (option == 8):
			plobj.format = newValue.strip(" ")
			if (not (plobj.format in availableFormats)):
				sys.stderr.write("\033[91mAvailable formats are %s.\n\033[0m" % ", ".join(availableFormats))
		else:
			print("\033[91mInvalid option.\033[0m")

	return None


def PromptPlotParameters(plobj):
	# Prompt the used for plot parameters or modify the currently available plot parameters
	print("\033[95mFollowing are the current plot settings, select (option number;new parameter value). Press enter to finalize.\033[0m")
	print("\033[93m0. Don't modify anything\033[0m")
	print("\033[93m1. Nucleotide lengths to be plotted: \033[92m%s.\033[0m" % ", ".join(map(str, plobj.nucLengths)))
	print("\033[93m2. Plot title: \033[92m%s.\033[0m" % plobj.title)
	print("\033[93m3. X-label: \033[92m%s.\033[0m" % plobj.xLabel)
	print("\033[93m4. Y-label: \033[92m%s.\033[0m" % plobj.yLabel)
	print("\033[93m5. X-axis scale: \033[92mfrom %g to %g, in steps of %g\033[0m" % (plobj.xscale[0], plobj.xscale[2], plobj.xscale[1]))
	print("\033[93m6. Y-axis scale: \033[92mfrom %g to %g, in steps of %g\033[0m" % (plobj.yscale[0, 0], plobj.yscale[0, 2], plobj.yscale[0, 1]))
	print("\033[93m7. Advanced options\033[0m")

	completed = 0
	while completed == 0:
		print(">>"),
		stdout.flush()
		choices = input().strip("\n").strip(" ")
		if len(choices) > 1:
			(optionStr, newValue) = choices.split(";")
			option = int(optionStr.strip(" "))
		else:
			exitOption = int(choices.strip(" "))
			if (exitOption == 7):
				option = 7
			else:
				print("\033[2mParameters remain in thier current state.\033[0m")
				option = 0

		if (option == 0):
			completed = 1
		elif (option == 1):
			# Find the indices of the selected lengths in the forward and reverese matching data
			newLengths = map(lambda li: int(li.strip(" ")), newValue.split(","))
			plobj.lengths[0] = np.nonzero(np.in1d(plobj.forwardMatchData[:, 0], newLengths))
			plobj.lengths[1] = np.nonzero(np.in1d(plobj.reverseMatchData[:, 0], newLengths))

		elif (option == 2):
			plobj.title = newValue.strip(" ")
		elif (option == 3):
			plobj.xLabel = newValue.strip(" ")
		elif (option == 4):
			plobj.yLabel = newValue.strip(" ")
		elif (option == 5):
			plobj.xscale = map(lambda xi: int(xi.strip(" ")), newValue.split(","))
			# The range is exclusive of the right end point. So we must advance the end point by the increment step to include it in the range.
			pltObj.xscale[1] = pltObj.xscale[1] + pltObj.xscale[2]
			# Consistency checks
			nPts = np.arange(*plobj.xscale).shape[0]
			if (nPts < 2):
				sys.stderr.write("\033[2mThe scale will only contain %d points.\033[0m" % (nPts))
		elif (option == 6):
			newScale = map(lambda yi: float(yi.strip(" ")), newValue.split(","))

			# If the y scale is negative, we will interprett it as the scale for reverse mapping and if positive, we will interpret it as the scale for forward matching.
			# If it has both negative and positive numbers, we will interpret it as an axes scale for the combined matching
			for di in range(1):
				# If the scale is on the positive y-axis, it corresponds to forward matching
				# If the scale is for the negative y-axis, it corresponds to reverse matching
				if ((newScale[0] > (-1)**di) and (newScale[1] > (-1)**di)):
					plobj.yscale[di, :] = newScale


			nPts = (newScale[1] - newScale[0])/float(newScale[2])
			if (np.any(nPts < 2)):
				sys.stderr.write("\033[2mOne of the scales contain only %d points.\n\033[0m" % (nPts))
			if (np.array_equal(np.sign(plobj.yscale[:, 0]), np.sign(plobj.yscale[:, -1])) == False):
				sys.stderr.write("\033[2mY-axis scale seems to include 0, this could show mRNA indices that were not involved in the binding.\n\033[0m")

		####################################
		# Advanced options
		elif (option == 7):
			PromptAdvancedParameters(plobj)
		else:
			print("\033[91mInvalid option.\033[0m")

	return None


def Plot(plobj, isScaled = 0, quiet = 0):
	# Produce plots with the selected data files and the plot parameters
	print("Forward\n{}\nReverse\n{}".format(plobj.nForw, plobj.nRev))
	xreps = 1
	shift = 0
	plotfname = ("./../plots/%s_%s.pdf" % (plobj.pstvdFname[:plobj.pstvdFname.index(".")], plobj.sRNAPoolFname[:plobj.sRNAPoolFname.index(".")]))
	xaxis = np.arange(plobj.pstvdLength + shift)
	with PdfPages(plotfname) as pdf:
		for l in range(plobj.lengths.shape[0]):
			for d in range(2):
				if (d == 0):
					nuclens = np.where(np.in1d(plobj.nForw[:, 0], plobj.lengths[l]))[0]
					yaxis = np.tile(np.sum(plobj.nForw[nuclens, 1:], axis = 0), xreps)[:len(xaxis)]
				else:
					nuclens = np.where(np.in1d(plobj.nRev[:, 0], plobj.lengths[l]))[0]
					yaxis = (-1) * np.tile(np.sum(plobj.nRev[nuclens, 1:], axis = 0), xreps)[:len(xaxis)]
				
				# If there is a matching at position i, for length l nucleotide, then we want to set: y[i + j] = max(y[i + j], y[i]).
				print("yaxis\n{}\nlength\n{}".format(yaxis, plobj.lengths[0][l]))
				scaled_yaxis = np.zeros_like(yaxis)
				for i in range(yaxis.shape[0] - plobj.lengths[0][l]):
					if (yaxis[i] > 0):
						for j in range(plobj.lengths[0][l]):
							scaled_yaxis[i + j] = max(yaxis[i], yaxis[i + j])
				print("scaled yaxis\n{}\nlength\n{}".format(scaled_yaxis, plobj.lengths[0][l]))
				# yaxis = scaled_yaxis

				fig = plt.figure(figsize = (26, 18))
				# plt.plot(xaxis, yaxis, color = "0.5")
				currax = plt.gca()
				currax.set_xlabel(plobj.pstvdFname, fontsize = 36)
				currax.set_ylabel(plobj.sRNAPoolFname, fontsize = 36)
				rightax = currax.twinx()
				topax = currax.twiny()
				if (d == 0):
					currax.set_ylim([0, 5])
					rightax.set_ylim([0, 5])
					currax.plot(xaxis, scaled_yaxis, marker="o", markersize=30, alpha=0.75, color="red", linestyle="None")
					# Draw a line from at y, from x to x + l
					# for i in range(len(xaxis)):
					# 	if (yaxis[i] > 0):
					# 		currax.plot([xaxis[i], xaxis[i] + plobj.lengths[0][l] - 1], [yaxis[i], yaxis[i]], color="red", linestyle="-", linewidth=5)
				else:
					currax.set_ylim([-5, 0])
					rightax.set_ylim([-5, 0])
					topax.plot(xaxis, yaxis, marker="o", markersize=30, alpha=0.75, color="red", linestyle="None")
					# for i in range(len(xaxis)):
					# 	if (yaxis[i] > 0):
					# 		topax.plot([xaxis[i], xaxis[i] + plobj.lengths[0][l] - 1], [yaxis[i], yaxis[i]], color="red", linestyle="-", linewidth=5)
				currax.set_xticks(xaxis, [tc for tc in range(1, 1 + len(xaxis))])
				topax.set_xticks(xaxis, [tc for tc in range(1, 1 + len(xaxis))])
				# currax.set_yticks(np.linspace(0, 20000, 10))
				# currax.grid(which = 'both', color = '0.9', linestyle = '-', linewidth = 2)

				plt.title("direction = %d, length %s" % (d, np.array_str(plobj.lengths[l])), fontsize = 48, y = 1.05)
				print("xaxis\n{}\nnonzero\n{}".format(xaxis, np.nonzero(scaled_yaxis)))
				
				plt.fill_between(xaxis, np.zeros(scaled_yaxis.shape[0], dtype = np.float), scaled_yaxis, color = "0.5")
				
				currax.tick_params(axis='y', which='both', pad = 20, direction = 'out', length = 10, width = 3)
				# rightax.tick_params(axis='y', which='both', pad = 20, left = 'off', right = 'off', labeltop = 'off', labelbottom = 'off', labelright = 'off', labelleft = 'off', direction = 'in', length = 10, width = 3)
				if (d == 0):
					currax.tick_params(axis='x', which='both', pad = 20, direction = 'out', length = 10, width = 3)
					topax.tick_params(axis='x', which='both', pad = 20, direction = 'in', length = 10, width = 3)
				else:
					currax.tick_params(axis='x', which='both', pad = 20, direction = 'in', length = 10, width = 3)
					topax.tick_params(axis='x', which='both', pad = 20, direction = 'out', length = 10, width = 3)
				# topax.set_xticks(currax.get_xticks())
				# rightax.set_yticks(currax.get_yticks())
				
				print("Saving plot for length %s." % (np.array_str(plobj.lengths[l])))
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


def Normalize(plotObjs, genes, pools, tols):
	# Multiply the number of reads at every position by a constant that depends on the nucleatide length.
	# This is to ensure that all reads are scaled for a pool size that is 1,000,000.
	# Scaling array provides the scaling constants by pool gene tol length.
	# refer by: scaling[pool][gene][tol][nuclen]
	scaling = {"pool":{"gene":{0:{0:0}}}}
	with open("scalings.txt", 'r') as sfp:
		for line in sfp:
			if (not (line[0] == "#")):
				contents = map(lambda cont: cont.strip(" "), line.strip("\n").strip(" ").split(" "))
				# print "contents"
				# print contents
				if (not (contents[0] in scaling)):
					scaling.update({contents[0]:{}})
				if (not (contents[1] in scaling[contents[0]])):
					scaling[contents[0]].update({contents[1]:{}})
				if (not (contents[2] in scaling[contents[0]][contents[1]])):
					scaling[contents[0]][contents[1]].update({contents[2]:{}})
				if (not (contents[3] in scaling[contents[0]][contents[1]][contents[2]])):
					scaling[contents[0]][contents[1]][contents[2]].update({int(contents[3]):[np.float(contents[4]), np.float(contents[5])]})

	# print "pool"
	# print pools
	# print "genes"
	# print genes
	# print "tols"
	# print tols
	print("scaling")
	print(scaling)

	# Write values into the normalized datasets
	for i in range(len(pools)):
		if (pools[i] in scaling.keys()):
			for j in range(len(genes)):
				if (genes[j] in scaling[pools[i]].keys()):
					for t in range(len(tols)):
						if (tols[t] in map(int, scaling[pools[i]][genes[j]].keys())):
							# print("lengths\n%s" % (np.array_str(plotObjs[i][j].forwardMatchData[:, 0])))
							# print("scaling[pool][genes][tols]")
							# print scaling[pools[i]][genes[j]][str(tols[t])]
							for l in range(plotObjs[i][j].forwardMatchData.shape[0]):
								if (plotObjs[i][j].forwardMatchData[l, 0] in scaling[pools[i]][genes[j]][str(tols[t])].keys()):
									# print "plotObjs[i][j].forwardMatchData[l, 1:]"
									# print plotObjs[i][j].forwardMatchData[l, 1:]
									# print "scaling[pools[i]][genes[j]][str(tols[t])][plotObjs[i][j].forwardMatchData[l, 0]]"
									# print scaling[pools[i]][genes[j]][str(tols[t])][plotObjs[i][j].forwardMatchData[l, 0]]
									plotObjs[i][j].nForw[l, 0] = plotObjs[i][j].forwardMatchData[l, 0]
									plotObjs[i][j].nForw[l, 1:] = plotObjs[i][j].forwardMatchData[l, 1:] * scaling[pools[i]][genes[j]][str(tols[t])][plotObjs[i][j].forwardMatchData[l, 0]][0]
									plotObjs[i][j].nRev[l, 0] = plotObjs[i][j].reverseMatchData[l, 0]
									plotObjs[i][j].nRev[l, 1:] = plotObjs[i][j].reverseMatchData[l, 1:] * scaling[pools[i]][genes[j]][str(tols[t])][plotObjs[i][j].reverseMatchData[l, 0]][1]
									# print("%s" % (np.array_str(plotObjs[i][j].nForw[l, :], max_line_width = 150)))
					# Write the normalized matching data into file
					WriteMatrixToFile(plotObjs[i][j].scaledForwDataFile, plotObjs[i][j].nForw, 0, 0, 0, dataType = 'f')
					WriteMatrixToFile(plotObjs[i][j].scaledRevDataFile, plotObjs[i][j].nRev, 0, 0, 0, dataType = 'f')
	return None


def ReadListFromFile(fname, ignore = 0):
	# Read a list from a file when the number of rows and columns is not directly specified
	# the first few lines can be ignored. The number of such lines is specified as "ignore"
	data = []
	ignored = []
	with open(fname, 'r') as fid:
		for (li, line) in enumerate(fid):
			contents = map(int, line.strip("\n").strip(" ").split(" "))
			if (li >= ignore):
				data.append(contents)
			else:
				ignored.append(contents)
	
	# convert the data into an np array
	structred = np.array(data, dtype = int)
	return (structred, ignored)


def CompileLaTex(fname):
	# Compile a tex file using pdflatex
	texCompileCommand = "/Library/TeX/texbin/pdflatex " + fname + " >/dev/null"
	os.system(texCompileCommand)
	# Remove the aux and log files
	os.remove(fname[:-3] + "aux")
	os.remove(fname[:-3] + "log")
	return None


def Report(plobj):
	# Write the list of all sRNA sequences that matched at a given index on the mRNA sequence
	# Read the raw forward and reverse matching files associated to this matching and write the list of sRNA sequences whose index is recorded.
	rawForwMatchingFname = ("raw_forward_matchings_%s_%s_tol%d.txt" % (plobj.sRNAPoolFname[:-4], plobj.pstvdFname[:-4], plobj.tolerance))
	rawRevMatchingFname = ("raw_reverse_matchings_%s_%s_tol%d.txt" % (plobj.sRNAPoolFname[:-4], plobj.pstvdFname[:-4], plobj.tolerance))

	(rawForwData, __)  = ReadListFromFile(rawForwMatchingFname, ignore = 1)
	# (rawRevData, __)  = ReadListFromFile(rawRevMatchingFname, ignore = 1)

	# Sort the data according to the second column -- which contains the indices of nucleotides which are involved in matchings
	rawForwData = rawForwData[np.argsort(rawForwData[:, 0], kind = "quicksort")]
	nForwMatches = rawForwData.shape[0]
	# rawRevData = rawRevData[np.argsort(rawRevData[:, 0], kind = "quicksort")]
	# nRevMatches = rawRevData.shape[0]

	# Read the sRNA sequences that correspond to the indices recorded in the raw-matching arrays
	ExplicitSequences = {mRNAIdx:{} for mRNAIdx in np.unique(rawForwData[:, 1])}
	nNucs = 0
	forwMatchCount = 0
	lineMapping = []
	with open(plobj.sRNAPoolFname, 'r') as sRNAFid:
		for (lNo, line) in enumerate(sRNAFid):
			if ((lNo % 4) == 1):
				sRNASeq = line.strip("\n").strip(" ")		
				if (not ("N" in sRNASeq)):
					lineMapping.append([nNucs, lNo])
					if (forwMatchCount < nForwMatches):
						while (nNucs == rawForwData[forwMatchCount, 0]):
							sRNALen = len(sRNASeq)
							print("Sequence %s of length %d at position %d of the raw matching info has a match with gene position %d." % (sRNASeq, sRNALen, forwMatchCount, rawForwData[forwMatchCount, 1]))
							if (not (sRNALen in ExplicitSequences[rawForwData[forwMatchCount, 1]])):
								ExplicitSequences[rawForwData[forwMatchCount, 1]].update({sRNALen:[]})
							ExplicitSequences[rawForwData[forwMatchCount, 1]][sRNALen].append(sRNASeq)
							forwMatchCount = forwMatchCount + 1
							if (forwMatchCount == nForwMatches):
								break
					nNucs = nNucs + 1

	# Debugging purposes
	WriteMatrixToFile("mapping.txt", np.array(lineMapping, dtype = int), 0);

	# Write the Explicit sequences dictionary to a file, for future plotting
	
	# Write the list of sequences to a latex file in the following manner.
	# The first page will have the full mRNA sequence in with a number to each index. Each nucleotide in the mRNA sequence will be a hyperlink to the page that contains the list of all sRNA sequences (lengthwise) that match with the particular mRNA sequence.
	# Read the mRNA sequence
	with open(plobj.pstvdFname, 'r') as mRNAFid:
		mRNASequence = mRNAFid.readline().strip("\n").strip(" ")

	explicitMatchFname = ("explicit_forward_%s_%s_tol%d.tex" % (plobj.sRNAPoolFname[:-4], plobj.pstvdFname[:-4], plobj.tolerance))

	with open(explicitMatchFname, 'w') as expFid:
		# Latex preamble
		expFid.write("\\documentclass[a4paper, 12pt]{article}\n")
		expFid.write("\\usepackage[left = 2cm, right = 2cm, top = 2cm, bottom = 2cm, landscape]{geometry}\n")
		expFid.write("\\usepackage{float}\n")
		expFid.write("\\usepackage{getfiledate,booktabs}\n")
		expFid.write("\\usepackage[pdftex,pdfauthor={Pavithran S Iyer}, pdftitle={\\bf RNA Silencing}]{hyperref}\n")
		expFid.write("\hypersetup{colorlinks,citecolor=blue,filecolor=blue,linkcolor=blue,urlcolor=blue}\n")

		expFid.write("\\title{Explicit sRNA sequences that bind to a particular gene on the mRNA sequence}\n")
		expFid.write("\\author{Pavithran Iyer}\n")
		expFid.write("\\begin{document}\n")
		expFid.write("\\maketitle\n")
		
		# Write the details of the time and date where the report was created
		expFid.write("\\section*{File details}\n")
		expFid.write("{\Large\sf\n")
		expFid.write("\\getfiledate[head=\\baselineskip,prefix=The date of final changes to file,foot=1ex,marker={$\\blacktriangleright$},markercolor=purple,filenamecolor=blue,width=0.95\\textwidth,datecolor=red,inlinespace=.5em,align=left,boxed,separator=at,framesep=1cm,framerule=2pt,framecolor=olive!25]")
		
		fnameTeXFormat = explicitMatchFname.replace("_", "\\string_")
		expFid.write("{%s}\n" % fnameTeXFormat)

		expFid.write("\\vspace{1cm}\n")

		expFid.write("Matching details:\n")
		expFid.write("\\begin{table}[H]\n")
		expFid.write("\\Large\\tt\n")
		expFid.write("\\begin{center}\n")
		expFid.write("\\begin{tabular}{c|c}\n")
		expFid.write("\\toprule\n")
		expFid.write("mRNA sequence & %s \\\\ \\midrule\n" % (plobj.pstvdFname.replace("_", "\\_")))
		expFid.write("sRNA Pool & %s \\\\ \\midrule\n" % plobj.sRNAPoolFname.replace("_", "\\_"))
		expFid.write("Number of mismatches allowed & %d \\\\ \\bottomrule\n" % (plobj.tolerance))
		expFid.write("\end{tabular}\n")
		# expFid.write("\caption{Details of the simulation.}\n")
		expFid.write("\end{center}\n")
		expFid.write("\end{table}\n")
		
		expFid.write("\n\\clearpage\n\n")

		for (mi, nuc) in enumerate(mRNASequence):
			if (mi in ExplicitSequences):
				expFid.write("\\parbox{0.4cm}{%d \\\\ " % mi)
				expFid.write("\\hyperref[mRNA%d]{%s}}" % (mi, nuc))
			else:
				expFid.write("\\parbox{1cm}{%d \\\\ " % mi)
				expFid.write("%s}" % (nuc))

			expFid.write("\n")

		for mrnai in ExplicitSequences:
			expFid.write("\\section{Sequences that matched at index %d} \\label{mRNA%d}\n" % (mrnai, mrnai))
			if (ExplicitSequences[mrnai]):
				for li in ExplicitSequences[mrnai]:
					expFid.write("\\subsection{sRNA sequences with %d nucleotides}\n" % li)
					for si in range(len(ExplicitSequences[mrnai][li])):
						expFid.write("%s\n" % (ExplicitSequences[mrnai][li][si]))
				expFid.write("\n\\clearpage\n\n")

		expFid.write("\\end{document}")
	return None


def Juxtapose(referencePlot, plotToBeScaled, direction):
	# Compare two plots by scaling the data of one plot to fit another plot
	print("\033[95mAccording to which crieteria do you want to juxtapose plots %s and %s?\033[0m" % (referencePlot.name, plotToBeScaled.name))
	print("\033[93m1 -- Equal maximum\n2 -- Equal total reads.\033[0m")
	print(">>"),
	crieteria = int(input().strip("\n"))

	if crieteria == 1:
		if direction == 0:
			basePlotMax = max(referencePlot.forwardMatchData[:, 0])
			scalingMax = max(plotToBeScaled.forwardMatchData[:, 0])
		else:
			basePlotMax = max(referencePlot.reverseMatchData[:, 0])
			scalingMax = max(plotToBeScaled.reverseMatchData[:, 0])

		scale = basePlotMax/float(scalingMax)
		print("\033[2mAll reads in %s will be multiplied by a factor of %.2e.\033[0m" % (plotToBeScaled.name, scale))

		if direction == 0:
			for pi in range(plotToBeScaled.forwardMatchData.shape[0]):
				plotToBeScaled.forwardMatchData[pi, 0] = plotToBeScaled.forwardMatchData[pi, 0] * scale
		else:
			for pi in range(plotToBeScaled.reverseMatchData.shape[0]):
				plotToBeScaled.reverseMatchData[pi, 0] = plotToBeScaled.reverseMatchData[pi, 0] * scale

	elif crieteria == 2:
		if direction == 0:
			totalReadsInRef = np.sum(referencePlot.forwardMatchData[:, 0])
			totalReadsInChange = np.sum(plotToBeScaled.forwardMatchData[:, 0])
			print("\033[2mEuqalizing total reads in %s to %d.\033[0m" % (plotToBeScaled.name, totalReadsInRef))
			for pi in range(plotToBeScaled.forwardMatchData.shape[0]):
				plotToBeScaled.forwardMatchData[pi, 0] = plotToBeScaled.forwardMatchData[pi, 0] * scale
		else:
			totalReadsInRef = np.sum(referencePlot.reverseMatchData[:, 0])
			for pi in range(plotToBeScaled.reverseMatchData.shape[0]):
				plotToBeScaled.reverseMatchData[pi, 0] = plotToBeScaled.reverseMatchData[pi, 0] * scale
	else:
		print("\033[95mCancelled.\033[0m")

	return None

if __name__ == '__main__':
	completed = 0
	isMenu = 1
	canvanses = {}
	# list of plot data for which normalized versions are available
	scaledAvailable = []
	while completed == 0:
		if (isMenu == 1):
			print("\033[92m**** MENU ****\033[0m")
			print("\033[93m0 -- Quit\033[0m")
			print("\033[93m1 -- Load new data for matching\033[0m")
			print("\033[93m2 -- Plot a previously loaded data\033[0m")
			print("\033[93m3 -- Relative scaling of matching data\033[0m")
			print("\033[93m4 -- Compare data and plot\033[0m")
			print("\033[93m5 -- Save matching data to a file\033[0m")
			print("\033[93m6 -- Report")
			print("\033[93m7 -- Menu")
			isMenu = 0
		else:
			print("\033[92mDone, what next? choose from Menu.\033[0m")

		print(">>"),
		userChoice = int(input().strip("\n").strip(" "))

		if (userChoice == 0):
			completed = 1
		elif (userChoice == 1):
			inputs = [("example_gene.txt", "example_pool.txt", 0)]
			name = "default"
			for (gene, pool, tol) in inputs:
				newCan = Canvas(gene, pool, tol)
				Load(newCan)
				canvanses.update({name:newCan})

		elif (userChoice == 2):
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
			userSelection = input().strip("\n").strip(" ")
			if (userSelection == "all"):
				plotDataChoice = range(len(canvanses.keys()))
			else:
				plotDataChoice = [int(userSelection)]
			
			for pdi in plotDataChoice:
				Plot(canvanses[list(canvanses.keys())[pdi]], isScaled = 1, quiet = 1)

		elif (userChoice == 3):
			print("\033[93mOf the %d available matching data sets, which ones do you want to include?\033[0m" % len(canvanses.keys()))
			for (ni, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (ni, name))
			print(">>"),
			userSelection = input().strip("\n").strip(" ")
			
			if (userSelection == "all"):
				normDataSets = range(len(canvanses.keys()))
			else:
				normDataSets = map(int, userSelection.split(","))

			normObjects = [canvanses[canvanses.keys()[ds]] for ds in normDataSets]

			genes = list(set([plo.pstvdFname[:plo.pstvdFname.index(".")] for plo in normObjects]))
			pools = list(set([plo.sRNAPoolFname[:plo.sRNAPoolFname.index(".")] for plo in normObjects]))
			tols = [0]
			
			normMatrix = [[None for pl in pools] for gn in genes]
			for li in range(len(normObjects)):
				gidx = genes.index(normObjects[li].pstvdFname[:normObjects[li].pstvdFname.index(".")])
				pidx = pools.index(normObjects[li].sRNAPoolFname[:normObjects[li].sRNAPoolFname.index(".")])
				normMatrix[gidx][pidx] = normObjects[li]

			# plotting is not included.
			Normalize(normMatrix, genes, pools, tols)
			
			normDataNames = [canvanses.keys()[ds] for ds in normDataSets]
			for di in range(len(normDataNames)):
				if not (normDataNames[di] in scaledAvailable):
					scaledAvailable.append(normDataNames[di])

		elif (userChoice == 4):
			print("\033[93mOf the %d available matching data sets, which ones would you like to compare and modify?\033[0m" % len(canvanses.keys()))
			for (ni, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (ni, name))

			print(">>"),
			compareChoices = map(int, input().strip("\n").strip(" ").split(",")[:2])

			plotsToCompare = [canvanses[canvanses.keys()[compareChoices[0]]], canvanses[canvanses.keys()[compareChoices[1]]]]

			isReverse = 1
			isForward = 1
			for pobj in plotsToCompare:
				if (np.sum(pobj.reverseMatchData) == 0):
					isReverse = 0
				if (np.sum(pobj.forwardMatchData) == 0):
					isForward = 0

			if (isReverse * isForward > 0):
				print("\033[93mCompare forward or reverse mapping data?\n0 -- Forward\n1 -- Reverse\033[0m")
				print(">>"),
				compareDir = int(input().strip("\n").strip(" "))
			elif (isReverse + isForward == 0):
				print("\033[91mMatching data is void. Cannot compare.\033[0m")
			elif isReverse == 0:
				compareDir = 0
			elif isForward == 0:
				compareDir = 1
			else:
				pass
			Juxtapose(plotsToCompare[0], plotsToCompare[1], compareDir)
		
		elif (userChoice == 5):
			print("\033[93mOf the %d available matching data sets, which one would you like to plot?\033[0m" % len(canvanses.keys()))
			for (ni, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (ni, name)),
				if (name in scaledAvailable):
					print("\033[92m\t (*)\033[0m")
				else:
					print("")

			print(">>"),
			saveDataChoice = int(input().strip("\n").strip(" "))

			print("\033[95m!! This will overwrite the existing data for %s.\n0 -- Cancel\n1 -- OK.\033[0m" % canvanses.keys()[saveDataChoice])
			print(">>"),
			isOverwrite = int(input().strip("\n").strip(" "))
			if (isOverwrite == 1):
				Save(canvanses[canvanses.keys()[saveDataChoice]])

		elif (userChoice == 6):
			print("\033[93mOf the %d available matching data sets, for which one would you like to a report?\033[0m" % len(canvanses.keys()))
			for (ni, name) in enumerate(canvanses):
				print("\033[93m%d -- %s\033[0m" % (ni, name))
			
			print(">>"),
			reportDataChoice = int(input().strip("\n").strip(" "))
			Report(canvanses[canvanses.keys()[reportDataChoice]])

		elif (userChoice == 7):
			isMenu = 1

		else:
			print("\033[91mUnknown option!\033[0m")

	print("\033[92mDone.\033[0m")
