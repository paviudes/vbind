import datetime as dt
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

class Data:
	"""
	Class describing the properties of the matching data that is required for plotting.
	"""
	def __init__(self, gene, pool, lengths, tol):
		self.lengths = np.array(lengths, dtype = np.int)
		self.gene = gene
		self.pool = pool
		self.tol = tol
		self.forward = np.empty([])
		self.reverse = np.empty([])

	def View(self):
		# View all the parameters of the class
		print("gene = {}".format(self.gene))
		print("pool = {}".format(self.pool))
		print("tol = {}".format(self.tol))
		print("lengths = {}".format(self.lengths))
		print("forward = {}".format(self.forward))
		print("reverse = {}".format(self.reverse))
		return None

	def Update(self, prop, value):
		# Update a particular property.
		exec("self.%s = %s" % (prop, value))
		return None

	def Derived(self):
		# Load the derived parameters.
		with open("./../data/input/%s.txt" % self.gene, "r") as gf:
			self.genelen = len(gf.readline().strip("\n").strip(" "))
		return None

	def ForwardMatchingFile(self):
		# Name of the file containing the forward matching reads
		fname = ("./../data/output/forward_matchings_%s.f_%s_tol%d.txt" % (self.pool, self.gene, self.tol))
		return fname

	def ReverseMatchingFile(self):
		# Name of the file containing the reverse matching reads
		fname = ("./../data/output/reverse_matchings_%s.f_%s_tol%d.txt" % (self.pool, self.gene, self.tol))
		return fname

	def Collect(self):
		# Collect the forward and reverse matchings data
		self.forward = self.ParseReads(self.ForwardMatchingFile())
		self.reverse = self.ParseReads(self.ReverseMatchingFile())
		return None

	def ParseReads(self, rfname):
		# Open the file containing the reads and add up the reads according to the prescription encoded in lengths.
		reads = np.zeros((self.lengths.shape[0], self.genelen), dtype = np.int)
		for i in range(reads.shape[0]):
			with open(rfname, "r") as rfp:
				for line in rfp:
					data = list(map(np.int, line.strip("\n").strip(" ").split(" ")))
					if (data[0] in self.lengths):
						reads[i, :] = reads[i, :] + data[1:]
		return reads

	def Help(self):
		# Mannual containing the class variables, their types and default values.
		self.mannual = mannual = {"lengths":["Neucleotide lengths for which the reads must be retrieved and plotted", "A list of lists. Every set of neucleotide lengths for which the reads must be added up and the sum plotted, must be given as a list. If several plots are required, one for each set of nuleotide lengths, those sets must be provided as a list of lists", "[[<minimum neucleotide length>, <average nucleotide length>, <maximum nucleotide length>]]"],
			   					  "gene":["Filename containing the gene. We don't need the file to even exist, just its name is used to contruct the name of the file containing the reads", "String filename with extension", "noname"],
			   					  "pool":["Filename containing the neucleotide pool. We don't need the file to even exist, just its name is used to contruct the name of the file containing the reads", "String filename with extension", "noname"],
			   					  "tol":["Number of mismatches", "Whole number > 0 indicating the number of mismatches with which the reads were collected", "0"]}
		for (kw, val) in enumerate(mannual):
			print("\"%s\"\n\tDescription: %s\n\tUsage:\n\t%s\n\tDefault:\n\t%s" % (kw, val[0], val[1], val[2]))
		return None


class Plot():
	"""
	Variables and functions concerning plotting the mathcing reads.
	"""
	def __init__(self, data):
		# Data associated to the plot
		self.data = data
		# Global options
		self.xreps = 1
		self.xshift = 0
		self.isnorm = 0
		self.iscircular = 0
		self.plotfname = "./../plots/%s_%s_t%d_%s.pdf" % (self.data.gene, self.data.pool, self.data.tol, "".join(list(map(lambda num: "%d" % num, self.data.lengths.ravel()))))
		# Limits
		self.xlims = [self.xshift, self.xreps * self.data.genelen]
		self.ylimsforward = [0, 1.1 * np.max(self.data.forward)]
		self.ylimsreverse = [-1.1 * np.max(self.data.reverse), 0]
		# Plot data
		self.xaxis = np.tile(np.linspace(0, self.data.genelen - 1, self.data.genelen, dtype = np.int), self.xreps)[self.xlims[0]:self.xlims[1]]
		self.yaxisforward = self.data.forward
		self.yaxisreverse = (-1) * self.data.reverse
		# Ticks
		self.xticks = self.xaxis[np.linspace(0, self.xaxis.shape[0] - 1, 10, dtype = np.int)]
		self.xticklabels = self.xaxis[np.linspace(0, self.xaxis.shape[0] - 1, 10, dtype = np.int)]
		self.yticksforward = np.linspace(self.ylimsforward[0], self.ylimsforward[1], 10, dtype = np.int)
		self.yticklabelsforward = np.linspace(self.ylimsforward[0], self.ylimsforward[1], 10, dtype = np.int)
		self.yticksreverse = np.linspace(self.ylimsreverse[1], self.ylimsreverse[0], 10, dtype = np.int)
		self.yticklabelsreverse = np.linspace(self.ylimsreverse[1], self.ylimsreverse[0], 10, dtype = np.int)
		self.tclength = 1
		self.tcwidth = 1
		# Colors and fills
		self.linecolor = 'k'
		self.topfill = "0.75"
		self.bottomfill = "0.25"
		self.linestyle = '-'
		# Axes labels
		self.xlabel = "Neucleotide index"
		self.ylabeltop = "Forward matching reads"
		self.ylabelbottom = "Forward matching reads"

	def View(self):
		# View the settings of plot.
		print(self.__dict__)

	def Update(self, prop, value):
		# Update a particular property.
		exec("self.%s = %s" % (prop, value))
		print("Property {} successfully modified to {}.".format(prop, value))
		return None

	def Help(self):
		# Display a cheat sheet of commands and their actions
		mannual = {"xlims":["Limits for X axis", "Either a range specification as: <low>,<high>.", "0,<max gene length>"],
				   "ylims":["Limits for X axis", "Either a range specification as: <low>,<high>.", "0,<max number of matches>"],
				   "linecolor":["Color of the plot curve", "Any color from, cf. https://matplotlib.org/2.0.2/api/colors_api.html", "0.75"],
				   "topfill":["Color to fill the plot area for forward matchings", "Any color from, cf. https://matplotlib.org/2.0.2/api/colors_api.html", "0.75"],
				   "bottomfill":["Color to fill the plot area for reverse matchings", "Any color from, cf. https://matplotlib.org/2.0.2/api/colors_api.html", "0.25"],
				   "xreps":["Number of times the X-axis needs to be repeated", "Whole number > 0", "1"],
				   "xshift":["Offset for the X-axis. In this case, the offset region will be plotted at the end if \"iscircular\" is turned on", "Whole number > 0", "0"],
				   "xticks":["Ticks on the X axis", "Either a range specification as: <low>,<high>,<number of points> or an explicit list of numbers.", "0,<upper limit in \"xlims\">","10"],
				   "yticks":["Ticks on the Y axis", "Either a range specification as: <low>,<high>,<number of points> or an explicit list of numbers.", "0,<upper limit in \"ylims\">","10"],
				   "tcwidth":["Thickness of the X and Y axis ticks, in points", "A number > 0", "1"],
				   "tclength":["Length of the X and Y axis ticks, in points", "A number > 0", "1"],
				   "xlabel":["Label on the X axis", "String", "Nucleotides"],
				   "ylabeltop":["Label for forward matchings on the Y-axis", "String", "Reads"],
				   "ylabelbottom":["Label for reverse matchings on the Y-axis", "String", "Reads"],
				   "issplit":["Should the plots for forward and reverse matchings must be split?", "1 for yes and 0 for no", "0"],
				   "isnorm":["Should the reads be normalized, so that the maximum number of reads is 1?", "1 for yes and 0 for no", "0"],
				   "iscircular":["Should the matchings be circular?", "1 for yes and 0 for no", "0"]}
		for (kw, val) in enumerate(mannual):
			print("\"%s\"\n\tDescription: %s\n\tUsage:\n\t%s\n\tDefault:\n\t%s" % (kw, val[0], val[1], val[2]))
		return None

	def RePlot(self):
		# Plot the data.
		with PdfPages(self.plotfname) as pdf:
			for l in range(self.data.lengths.shape[0]):
				for d in range(2):
					fig = plt.figure(figsize = (26, 18))
					ax = plt.gca()
					# For every length and direction (forward, reverse), create a new plot.
					if (d == 0):
						plt.plot(self.xaxis, self.yaxisforward[l], color = self.linecolor, linestyle = self.linestyle)
						plt.fill_between(self.xaxis, np.zeros_like(self.xaxis), self.yaxisforward[l], color = self.topfill)
						plt.ylim(self.ylimsforward)
						plt.yticks(self.yticksforward, self.yticklabelsforward)
					if (d == 1):
						plt.plot(self.xaxis, self.yaxisreverse[l], color = self.linecolor, linestyle = self.linestyle)
						plt.fill_between(self.xaxis, self.yaxisreverse[l], np.zeros_like(self.xaxis), color = self.bottomfill)
						plt.ylim(self.ylimsreverse)
						plt.yticks(self.yticksreverse, self.yticklabelsreverse)
					# X axis
					plt.xlim(self.xlims)
					plt.xticks(self.xticks, self.xticklabels)
					# font sizes of labels
					plt.tick_params(axis='both', which='both', labelsize = 36)
					# Saving the plot
					pdf.savefig(fig)
					plt.close()
			# Set the PDF attributes
			pdfInfo = pdf.infodict()
			pdfInfo['Title'] = ("Forward and Reverse matching reads")
			pdfInfo['Author'] = "Pavithran Iyer"
			pdfInfo['ModDate'] = dt.datetime.today()
			print("Plot modified and saved to %s." % (self.plotfname))
		return None


def LoadData(data, satisfied = 0):
	# Load the plot data.
	# Update the loaded data until user is satisfied.
	while (satisfied == 0):
		data.View()
		userinput = input(">>")
		if (userinput.lower == "done"):
			satisfied = 1
		else:
			(prop, value) = userinput.split(",")
			data.Update(prop, value)
	data.Derived()
	data.Collect()
	return None

def PlotReads(data, satisfied = 0):
	# Call the RePlot function until the user is satisfied of the plot.
	plot = Plot(data)
	plot.RePlot()
	while (satisfied == 0):
		# plot.View()
		userinput = input(">>").strip("\n")
		if (userinput.lower() == "done"):
			satisfied = 1
		else:
			(prop, value) = userinput.split(",")
			plot.Update(prop, value)
			plot.RePlot()
	return None

if __name__ == '__main__':
	# Load and plot reads data.
	pools = ["IIIweekCtrl_S4_R1_001-TR", "IIIweekCtrl_S4_R1_001-TR", "IIweekRDR_S3_R1_001-TR", "IIweekRDR_S3_R1_001-TR", "IIweekEV_S2_R1_001-TR", "IIweekEV_S2_R1_001-TR", "IIweekCtrl_S1_R1_001-TR", "IIweekCtrl_S1_R1_001-TR", "IIIweekRDR_S6_R1_001-TR", "IIIweekRDR_S6_R1_001-TR", "IIIweekEV_S5_R1_001-TR", "IIIweekEV_S5_R1_001-TR"]
	for pool in pools:
		for tol in range(2):
			data = Data("PSTVd_RG", pool, [[21, 22, 23, 24]], tol)
			LoadData(data, satisfied = 1)
			data.View()
			PlotReads(data, satisfied = 1)
