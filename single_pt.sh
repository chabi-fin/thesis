#!/bin/bash

mims=("FY1a" "FY1b" "FY2a" "FY2b" "FY3a" "FY3b" "FY5a" "FY5b" "FY6a" "FY6b")

for mim in "${mims[@]}"; do

    cd ~/thesis/stage2chg-param/
    cp run_g_job.sh extract_conformations.py ave_charges.py $mim
    cd ~/thesis/stage2chg-param/$mim
    
    mkdir -p configurations/
    cp extract_conformations.py configurations/
    cp fTyr_md.trr fTyr_md.gro run_g_job.sh configurations/
    cd configurations/
    python extract_conformations.py
    
    # make a list of the configuration subdirectories
    dirs=($(find . -type d -name "config_*"))
    home=$(pwd)
	
    # run single point calculation for every conformation 
    for dir in ${dirs[@]}
    do
	cp run_g_job.sh $dir
	cd $dir
	# remove all lines containing the string 'CONECT' from geo.pdb
	grep -F -v CONECT geo.pdb > geo.pdb.tmp && mv geo.pdb.tmp geo.pdb	
	
	antechamber -i geo.pdb -fi pdb -o geo.dat -fo gcrt -gv 1 -ge geo.gesp
	rm -f molecule.chk geo.gesp geo.esp	
	
	# edit calculation options
	sed -i -e 's/\opt//' geo.dat

	# set residue net charge
	if [[ $mim == "FY1a" || $mim == "FY2a" ]]; then
		sed -i -e 's/\0   1/-2  1/' geo.dat 
	else 
		sed -i -e 's/\0   1/-1  1/' geo.dat
	fi

	# run job for esp
	./run_g_job.sh
	
	cd $home

    done

    cd ~/thesis/stage2chg-param/

done

exit 0
