#!/bin/bash

MIMS=("FY1a" "FY1b" "FY2a" "FY2b" "FY3a" "FY3b" "FY5a" "FY5b" "FY6a")

for mim in "${MIMS[@]}"
do
	cd ~/thesis/
	cp -r amber14sb.ff stage2chg-param/$mim/
	cd stage2chg-param/
	cp *.mdp $mim/
	cp period $mim/
	cp residuetypes.dat $mim/
	cp run_job.sh $mim/
	cd $mim/
	

	# Generate topology using .pdb of capped residue
	gmxs pdb2gmx -ff amber14sb -f "$mim".pdb -o fTyr.gro -water tip4p -ignh -nobackup

	./run_job.sh &

done

exit 0
