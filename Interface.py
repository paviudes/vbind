import os
import sys
import time

class Submission():
	"""All the parameters required for a submission"""
	def __init__(self):
		self.sourceFiles = ["Binder.py", "run.sh", "post.sh"]
		self.bindingPairs = []
		self.tolerance = 0
		self.cores = 4
		self.job = "JOB"
		self.queue = "qwork"
		self.wall = 0

	def Display(self):
		# Display all the parameters and their current values
		print("\033[92m******************************************\033[0m")
		print("\033[93mjob --> \033[92m%s\033[0m" % (self.job))
		
		print("\033[93mBinding pairs --> \033[92m%d\033[0m" % (len(self.bindingPairs)))
		print("\033[92mgene\tpool\ttolerance\033[0m")
		for pi in range(len(self.bindingPairs)):
			print("\033[92m%s\t%s\t%d\033[0m" % (self.bindingPairs[pi][0], self.bindingPairs[pi][1], self.bindingPairs[pi][2]))

		print("")
		print("\033[93mcores --> \033[92m%d\033[0m" % (self.cores))
		print("\033[93mwall --> \033[92m%d\033[0m" % (self.wall))
		print("\033[93mqueue --> \033[92m%s\033[0m" % (self.queue))
		print("\033[92m******************************************\033[0m")

		quota = len(self.bindingPairs) * self.wall
		print("\033[2m%.3f%% of total usage quota will be used up if the simulation runs for the entrie walltime.\033[0m" % (quota/float(100000)))
		return None


def ConsoleInput():
	# Input from the console to set the parameters in the bqsubmit file
	submit = Submission()
	print("\033[93mEnter the pairs of genes and the sRNA pool which you want to bind.\nFormat must be <gene file name,pool file name,tolerance>\nWhen you are done, enter '"'Done'"'.\033[0m")
	done = 0
	nPairs = 0
	while (done == 0):
		print("\033[93m>>%d).\033[0m" % (nPairs + 1)),
		userInput = raw_input().strip("\n").strip(" ")
		if (userInput.lower() == "done"):
			done = 1
		else:
			(nucFile, poolFile, tol) = userInput.split(",")
			# Check if the files exist
			if (os.path.isfile(nucFile) and os.path.isfile(poolFile)):
				submit.bindingPairs.append([nucFile, poolFile, int(tol)])
				nPairs = nPairs + 1
			else:
				print("\033[91m\033[2mX Cannot include this combination, one the files are missing.\033[0m")

	print("\033[93mNumber of cores to be used at each node (default = 24).\033[0m")
	print("\033[93m>>\033[0m"),
	submit.cores = int(raw_input().strip("\n").strip(" "))

	print("\033[93mGive a name to the submission.\033[0m")
	print("\033[93m>>\033[0m"),
	submit.job = raw_input().strip("\n").strip(" ")

	print("\033[93mMaximum walltime per node in hours.\033[0m")
	print("\033[93m>>\033[0m"),
	submit.wall = int(raw_input().strip("\n").strip(" "))

	print("\033[93mName of the queue into which the submission must be launched.\nSelect one from '"'qwork'"', '"'qtest'"', '"'qfbb'"', '"'qfat256'"', '"'qfat512'"'.\nDefault is '"'qwork'"'.\033[0m")
	print("\033[93m>>\033[0m"),
	submit.queue = raw_input().strip("\n").strip(" ")

	return submit


def WriteToBQSubmit(submit):
	# Produce a bqsubmit.dat file from the parameters loaded into the submission object.
	# Produce the executable file for the simulation
	dateTimeCompleted = time.strftime("%d/%m/%Y %H:%M:%S")
	dateStamp = dateTimeCompleted.translate(None, "/: ")
	submit.sourceFiles = ["Binder.py", ("run_%s.sh" % dateStamp), ("post_%s.sh" % dateStamp)]
	with open(submit.sourceFiles[1], 'w') as runFid:
		# runFid.write("module add intel64/12.0.4.191 python64/2.7.5 multiprocess/0.70.4.dev0 mkl64/10.1.3.027\n")
		runFid.write("export OMP_NUM_THREADS=24\n")
		runFid.write("module load intel64/14.0.0.080 python64/2.7.6 mkl64/10.1.3.027\n")
		runFid.write("python Binder.py %d 2>&1 | tee -a log.txt\n" % submit.cores)
	os.system("chmod +x %s" % submit.sourceFiles[1])

	# Input parameter file
	with open("input.txt", 'w') as inpFid:
		inpFid.write("~~gene~~\n")
		inpFid.write("~~pool~~\n")
		inpFid.write("~~tol~~\n")

	# Post batch file
	resultsDir = ("%s_%s" % (submit.job, dateStamp))
	with open(submit.sourceFiles[2], 'w') as postFid:
		postFid.write("cat %s_*.BQ/input.txt >> SimulationInputs_%s.txt\n" % (submit.job, dateStamp))
		postFid.write("cat %s_*.BQ/details.txt >> SimulationResults_%s.txt\n" % (submit.job, dateStamp))
		
		postFid.write("mkdir %s\n" % resultsDir)
		postFid.write("cp -r SimulationResults_%s.txt SimulationInputs_%s.txt bqsubmit.dat %s/\n" % (dateStamp, dateStamp, resultsDir))
		
		# Raw matching files
		postFid.write("cp -r %s_*.BQ/raw_*.txt %s/\n" % (submit.job, resultsDir))
		# Forward matching files
		postFid.write("cp -r %s_*.BQ/forward_*.txt %s/\n" % (submit.job, resultsDir))
		# Reverse maching files
		postFid.write("cp -r %s_*.BQ/reverse_*.txt %s/\n" % (submit.job, resultsDir))

		postFid.write("tar -zcvf %s.tar.gz %s/\n" % (resultsDir, resultsDir))
	
	os.system("chmod +x %s" % submit.sourceFiles[2])

	geneFiles = [submit.bindingPairs[pi][0] for pi in range(len(submit.bindingPairs))]
	sRNAFiles = [submit.bindingPairs[pi][1] for pi in range(len(submit.bindingPairs))]
	tols = [submit.bindingPairs[pi][2] for pi in range(len(submit.bindingPairs))]

	linkFiles = geneFiles + sRNAFiles
	bindingParams = [("(%s, %s, %d)" % (submit.bindingPairs[pi][0], submit.bindingPairs[pi][1], submit.bindingPairs[pi][2])) for pi in range(len(submit.bindingPairs))]

	with open('bqsubmit.dat', 'w') as bqFid:
		bqFid.write("batchName = %s\n\n" % (submit.job))
		# link to the sRNA pool files
		bqFid.write("linkFiles = %s\n\n" % ", ".join(linkFiles))
		# copy the lattice text files
		bqFid.write("copyFiles = %s\n\n" % ", ".join(submit.sourceFiles))
		# Input file with the values of variables
		bqFid.write("templateFiles = input.txt\n\n")
		# Command to run on the compute node
		bqFid.write("command = ./%s\n\n" % (submit.sourceFiles[1]))
		# Number of jobs to be run simultaneously on a node
		bqFid.write("runJobsPerNode = 1\n\n")
		# Number of jobs to be accumulated on a node
		bqFid.write("accJobsPerNode = 1\n\n")
		# Command to be run after all jobs are exited or deleted
		bqFid.write("postBatch = ./post_%s.sh\n\n" % (dateStamp))
		# Required resources for each node
		bqFid.write("submitOptions=-q %s@mp2 -l walltime=%d:00:00,nodes=1\n\n" % (submit.queue, submit.wall))
		# Parameters to be scanned in the batch
		bqFid.write("param1 = (gene, pool, tol) = [%s]\n\n" % (", ".join(bindingParams)))
		# Number of nodes to be used concurrently
		bqFid.write("concurrentJobs = 144\n\n")
		# send an email
		bqFid.write("emailAddress = charith.adkar@gmail.com")
	return None

if __name__ == '__main__':
	submitObj = ConsoleInput()
	WriteToBQSubmit(submitObj)
	submitObj.Display()
