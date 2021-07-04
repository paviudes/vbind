import os
import numpy as np


def ReadMatrixFromFile(fname, dataType = 'i', atol = 1E-16):
	# Read a matrix from file
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


def ArrayToString(arr, data_type = 'i'):
	# Convert an array to string.
	if (data_type == 'i'):
		arr_str = ",".join(list(map(lambda x: "%d" % x, arr)))
	return arr_str

def IntArrayToString(arr):
	# Convert an integer array to string.
	return ",".join(list(map(lambda x: "%d" % x, arr)))

def ReverseComplement(seq):
	# Compute the reverse complement encoding.
	reverse_encoding = {"A":"T", "T":"A", "G":"C", "C":"G"}
	revcomp = [seq[s] for s in range(len(seq))]
	for s in range(len(seq)):
		revcomp[s] = reverse_encoding[seq[s]]
	return "".join(revcomp[::-1])

def LeastGreatestMultiple(number, factor):
	# Given a number, n, compute the lowest multiple, c, of a factor, x, such that: n <= c.x.
	# The multiple c = ceil(n/x)
	if (np.abs(factor) >= np.abs(number)):
		if number >= 0:
			return np.ceil(number)
		return np.floor(number)
	if number < 0:
		multiple = np.floor(number/factor)
		lcm = np.floor(multiple * factor)
	else:
		multiple = np.ceil(number/factor)
		lcm = np.ceil(multiple * factor)
	return lcm