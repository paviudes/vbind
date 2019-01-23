import os
# Read the input parameters in the new input format and for each input set create a text file containing its specifications.
with open("sim_data.txt", "r") as sf:
	for line in sf:
		(gene, pool, tol, cores) = line.strip("\n").strip(" ").split(" ")
		fname = "%s_%s_%s.txt" % (os.path.splitext(gene)[0], os.path.splitext(pool)[0], tol)
		with open(fname, "w") as of:
			of.write("%s\n%s\n%s" % (pool, gene, tol))