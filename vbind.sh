#!/bin/bash
# Check if the version of Python running is 3.x
version_info="$(python -c "import sys; print('{}.{}.{}'.format(*sys.version_info))")"
version_number="$(python -c "import sys; print('{}'.format(*sys.version_info))")"
default_check=1
python_command=python
if (( version_number <= 2 )); then
	python_command=python3
	default_check=0
fi
if [ "$default_check" -eq "0" ]; then
	if ! [ -x "$(command -v python3)" ]; then
		echo "Error: attempting to run using Python ${version_info}."
		echo "vbind needs Python 3 or higher to run."
		exit
	else
		echo -e "\033[2mSkipping the default Python ${version_info} installation and running ${python_command} instead.\033[0m"
	fi
fi

# Create the data output and plot folders.
mkdir -p data/output
mkdir -p data/plots

cd vbind
command=$(echo "${python_command} binder.py "$@"")
eval $command
cd ..