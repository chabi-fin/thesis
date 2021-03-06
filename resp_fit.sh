#!/bin/bash

# name and conformation of each mimetic
mims=("FY1a" "FY1b" "FY2a" "FY2b" "FY3a" "FY3b" "FY5a" "FY5b" "FY6a" "FY6b")

for mim in "${mims[@]}"; do

	# Set up list of directories to loop over
	cd ~/thesis/stage2chg-param/$mim/configurations/
	dirs=($(find . -type d -name "config_*"))
	home=$(pwd)
	
	for dir in ${dirs[@]}
	do
		# include input files
		cp resp1.in resp2.in resp1.qin $dir/
		cd $dir/	
		
		# process gaussian output to ESP
		if [ ! -f geo.gesp ]; then
		    echo $(pwd) 
		else
		    espgen -i geo.gesp -o geo.esp
		fi

		# first iteration of resp procedure
		resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg -q resp1.qin -e geo.esp

		# second iteration of resp procedure
		resp -O -i resp2.in -o resp2.out -p resp2.pch -t resp2.chg -q resp1.chg -e geo.esp
	
		cd $home
	done

done

exit 0
