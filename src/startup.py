def CheckDependencies():
	# Check if all the requires packages exist
	# if not, create a requirements text file.
	packages = [["scipy", "0.18", "Numerical operations", "Critical"],
				["numpy", "0.12", "Numerical operations", "Critical"],
				["multiprocessing", "1.0", "parallel computations", "Critical"],
				["matplotlib", "1.0", "plotting", "Critical"],
				["tqdm", "4.31.1", "progress bar", "Critical"],
				["datetime", "3.2", "Fetching date and time", "Mild"]]
	missing = []
	for i in range(len(packages)):
		try:
			exec("import %s" % (packages[i][0]))
		except:
			missing.append(packages[i])
	
	if (len(missing) > 0):
		print("\033[33m", end="")
		print("Missing or outdated packages might affect certain functionalities.")
		print("{:<10} | {:<10} | {:<30} | {:<10}".format("Package", "Version", "Affected functionality", "Impact"))
		print("{:<50}".format("---------------------------------------------------------------------"))
		
		is_critical = 0
		for i in range(len(missing)):
			if ("Critical" in missing[i][3]):
				is_critical = 1
				print("{:<10} | {:<10} | {:<30} | \033[7;31m{:<10}\033[0;33m".format(missing[i][0], missing[i][1], missing[i][2], missing[i][3]))
			else:
				print("{:<10} | {:<10} | {:<30} | \033[32m{:<10}\033[33m".format(missing[i][0], missing[i][1], missing[i][2], missing[i][3]))
		print("xxxxxx")
		print("\033[0m", end="")
		
		with open("./../requirements.txt", "w") as fp:
			fp.write("# Install the missing packages using pip install -r requirements.txt\n")
			for i in range(len(missing)):
				fp.write("%s>=%s\n" % (missing[i][0], missing[i][1]))
		print("\033[7;32mTo install all missing packages, run \"pip install -r requirements.txt\".\033[0m")

		if (is_critical):
			print("\033[0;31mExiting due to missing critical packages ...\033[0m")
			exit(0)
	return None


def DisplayLogoLicense():
	# Display logo as ascii drawing from http://ascii.mastervb.net with font = xcourb.tiff
	# Display license from the LICENSE file in chflow/
	logo = r"""
	       ###      ##            ###  
	        ##                     ##  
	### ##  ####   ###   ## ##   ####  
	 ## ##  ## ##   ##    ## ## ## ##  
	  ###   ## ##   ##    ## ## ## ##  
	  ###   ## ##   ##    ## ## ## ##  
	   #   #####  ######  ## ##  ##### 
	"""
	welcome = r"""
	Welcome to vbind version v1.0.
	Check out https://github.com/paviudes/vbind for help.
	"""
	license = r"""
	BSD 3-Clause License
	Copyright (c) 2021, Pavithran Iyer and Charith Adkar.
	All rights reserved.
	"""
	print("%s\033[0;36m%s%s\033[0m"% (logo, welcome, license))
	return None

def Usage():
	# print the usage of the binder.py script as well as the input file format.
	print("\033[2mUsage:\033[0m")
	print("\033[2mpython binder.py <input_file> <line_number>\033[0m")
	print("\033[2mwhere:\033[0m")
	print("\033[2m\t<input_file> should be a readable in data/input, formatted as follows.\033[0m")
	print("\033[2m\t\tEach line not beginning with \"#\" is: <gene> <pool> <tolerance> <circular> <cores>\033[0m")
	print("\033[2m\t\tsuch that\033[0m")
	print("\033[2m\t\t\t1. <gene>: name of the readable file in vbind/data/input that contains the gene sequence.\033[0m")
	print("\033[2m\t\t\t2. <pool>: name of the readable file in vbind/data/input that contains the pool.\033[0m")
	print("\033[2m\t\t\t3. <tolerance>: integer specifying the maximum number of mismatches allowed.\033[0m")
	print("\033[2m\t\t\t4. <circular>: binary (0 or 1) indicating circular (1) or linear (0) matching.\033[0m")
	# print("\033[2m\t\t\t5. <cores>: integer specifying the number of cores reserved for the script.\033[0m")	
	print("\033[2m\t\tAny line beginning with \"#\" is a comment and will be ignored by the script.\033[0m")
	print("\033[2m\t<line_number> should be an integer specifying the line number of the matching task in the input file.\033[0m")
	return None