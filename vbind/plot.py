import os
from sys import stdout
import datetime as dt
from tqdm import tqdm
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from utils import LeastGreatestMultiple, IntArrayToString

class PlotSettings:
	"""
	Settings used in plots.
	figsize: size of the plot frame
	fill: color for the fill or a grayscale number of the fill.
	font: font family for the text in the plot.
	fontsize: font size for the text in the plot.
	"""
	def __init__(self):
		self.figsize = (32, 24)

	def update(self, field, value):
		# Update the value of a field for plot settings.
		if (field == "figsize"):
			# Size of the frame for plot: must be a pair of integers separated by a comma.
			self.figsize = tuple(map(int, value.split(",")))
		elif (field == "fill"):
			# fill_shade: RG: 0.3, I: 0.65, example: 0.65
			# fill_shade = "0.3"
			self.fill = value
		elif (field == "font"):
			self.font = value
		elif (field == "fontsize"):
			self.fontsize = int(value)
		else:
			pass
		return None

	def load(self, fname):
		# Load plot settings from a file.
		with open(fname, "r") as fp:
			for line in fp:
				if (line[0] != "#"):
					(field, value) = line.strip("\n").strip(" ").split(" ")
					self.update(field.strip(" "), value.strip(" "))
		return None


def SetXTicks(xaxis):
	# Set the X-ticks.
	# If the length of the X-axis is less than 50, then we show an X-tick at every integer.
	# Else, we show an X-tick at every multiple of 50, besides the start and end positions.
	x_right = np.max(xaxis)
	if ((x_right > 0) and (x_right < 20)):
		xticks = np.arange(x_right + 1)
	elif ((x_right >= 20) and (x_right < 50)):
		xticks = np.unique(np.concatenate(([0], np.arange(4, x_right + 1, 5), [x_right + 1])))
	else:
		xticks = np.unique(np.concatenate(([0], np.arange(49, x_right + 1, 50), [x_right + 1])))
	xticklabels = list(map(lambda x: "%d" % (1 + x), xticks))
	return (xticks, xticklabels)


def SetYTicks(yaxis, direction):
	# Set the Y-ticks.
	# If the ticks is explicitly provided, simply use it, and set the Y-limit appropriately.
	# Else, do the following.
	# If the length of the Y-axis is:
	# 	1. less than 50, display a tick at every point.
	# 	2. between 50 and 500, display a tick at every 50.
	# 	3. between 500 and 1000, display a tick at every 100.
	# 	4. between 1000 and 5000, display a tick at every 500.
	# 	5. between 5000 and 15000, display a tick at every 1000.
	# 	6. greater than 15000, display a tick at every 5000.
	# In each of the cases, the Y-limit:
	# 	A. for forward matching is: LCM(1.1 * max(yaxis), tick interval)
	# 	B. for forward matching is: LCM(1.1 * min(yaxis), tick interval)
	yaxis_extent = np.max(np.abs(yaxis))
	if (yaxis_extent <= 50):
		if (direction == "forward"):
			yticks = np.arange(0, np.max(yaxis) + 1, 1)
		else:
			yticks = np.arange(0, np.min(yaxis) - 1, -1)
			# print("yaxis\n{}\nyticks: {}".format(yaxis, yticks))
	elif ((yaxis_extent > 50) and (yaxis_extent <= 500)):
		if (direction == "forward"):
			yticks = np.arange(0, np.max(yaxis) + 1, 50)
		else:
			yticks = np.arange(0, np.min(yaxis) - 1, -50)
	elif ((yaxis_extent > 500) and (yaxis_extent <= 1000)):
		if (direction == "forward"):
			yticks = np.arange(0, np.max(yaxis) + 1, 100)
		else:
			yticks = np.arange(0, np.min(yaxis) - 1, -100)
	elif ((yaxis_extent > 1000) and (yaxis_extent <= 5000)):
		if (direction == "forward"):
			yticks = np.arange(0, np.max(yaxis) + 1, 500)
		else:
			yticks = np.arange(0, np.min(yaxis) - 1, -500)
	elif ((yaxis_extent > 5000) and (yaxis_extent <= 15000)):
		if (direction == "forward"):
			yticks = np.arange(0, np.max(yaxis) + 1, 1000)
		else:
			yticks = np.arange(0, np.min(yaxis) - 1, -1000)
	else:
		# print("Set Yticks for length of the yaxis: {} and direction: {}".format(yaxis_extent, direction))
		if (direction == "forward"):
			yticks = np.arange(0, np.max(yaxis) + 1, round(np.max(yaxis)/10, -3))
		else:
			yticks = np.arange(0, np.min(yaxis) - 1, round(np.min(yaxis)/10, -3))

	# Tick labels
	yticklabels = list(map(lambda x: "%d" % (x), yticks))

	# Limits
	if (direction == "forward"):
		ylim = [0, LeastGreatestMultiple(1.05 * np.max(yaxis), 500)]
	else:
		ylim = [LeastGreatestMultiple(1.05 * np.min(yaxis), 500), 0]
		# print("yticks: {}\nylim: {}".format(yticks, ylim))

	return (yticks, yticklabels, ylim)


def Plot(plobj, plot_settings):
	# Produce plots with the selected data files and the plot parameters
	alignments = ["forward", "reverse"]
	plotfname = ("./../plots/%s_%s_tol_%d.pdf" % (plobj.genefname[:plobj.genefname.index(".")], plobj.poolfname[:plobj.poolfname.index(".")], plobj.tolerance))
	
	# Set the global parameters for the plot.
	matplotlib.rcParams["font.family"] = plot_settings.font

	with PdfPages(plotfname) as pdf:
		for l in range(len(plobj.lengths)):
			for d in range(2):
				fig = plt.figure(figsize = plot_settings.figsize)
				ax = plt.gca()
				if (d == 0):
					yaxis = plobj.scaled_forward[l, :]
				else:
					yaxis = plobj.scaled_reverse[l, :]
					ax.xaxis.tick_top()
				
				# Indices on the genome that have non-zero mappings.
				nonzero_gene_indicator = np.where(yaxis != 0, 1, 0)

				# plt.plot(xaxis, yaxis, color = "0.5")
				plt.fill_between(np.arange(plobj.gene_length), np.zeros(yaxis.shape[0], dtype = np.float), yaxis, where=nonzero_gene_indicator, color=plot_settings.fill)
				
				# Setting the X-ticks
				(xticks, xticklabels) = SetXTicks(np.arange(plobj.gene_length))
				ax.set_xticks(xticks)
				ax.set_xticklabels(xticklabels)

				# Setting the Y-ticks
				(yticks, yticklabels, ylim) = SetYTicks(yaxis, alignments[d])
				ax.set_yticks(yticks)
				ax.set_yticklabels(yticklabels)
				ax.set_ylim(ylim)

				plt.title("%s matching for lengths %s" % (alignments[d], IntArrayToString(plobj.lengths[l])), fontsize = 48, y = 1.1)
				
				ax.tick_params(axis='y', which='both', pad = 20, direction = 'inout', length = 10, width = 3, labelsize = plot_settings.fontsize)
				
				if (d == 0):
					ax.tick_params(axis='x', which='both', pad = 20, direction = 'inout', length = 10, width = 3, labelsize = plot_settings.fontsize)
				else:
					ax.tick_params(axis='x', which='both', pad = 20, direction = 'in', length = 10, width = 3, labelsize = plot_settings.fontsize)
				
				fig.tight_layout(pad=5)
				pdf.savefig(fig)  # saves the current figure into a pdf page
				plt.close()
			
			print("\033[2mDone plot for length %s.\033[0m" % (IntArrayToString(plobj.lengths[l])))
				
		# Set the PDF attributes
		pdfInfo = pdf.infodict()
		pdfInfo['Title'] = ("Comparing different sampling methods")
		pdfInfo['Author'] = "Pavithran Iyer"
		pdfInfo['ModDate'] = dt.datetime.today()
	print("\033[92mAll plots saved to %s.\033[0m" % plotfname)
	return None


def Visualize(rnap, lengths, is_scaled = 1, plot_settings = None):
	# Plot the results of the sRNA profiler.
	rnap.lengths = lengths
	GatherMatchingData(rnap, rnap.nForw, rnap.nRev)
	# Normalize the output if necessary
	if (is_scaled == 1):
		Normalize(rnap, rnap.gathered_forward, rnap.gathered_reverse)
		(rnap.gathered_forward, rnap.gathered_reverse) = (rnap.normalized_forward, rnap.normalized_reverse)
	ScaleMatchings(rnap, rnap.gathered_forward, rnap.gathered_reverse)
	# Load plot settings
	if plot_settings is None:
		plot_settings = PlotSettings()
		settings_fname = input(">>Settings file name (leave blank for default): ").strip("\n").strip(" ")
		if (len(settings_fname) == 0):
			settings_fname = "./../data/input/default_plot_settings.txt"
		plot_settings.load(settings_fname)
	# Plot
	Plot(rnap, plot_settings)
	return None