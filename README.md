# Thesis
Python programs and bash scripts related to my Masters thesis in AG Keller at FU Berlin

GET BONDED PARAMETERS FROM GAFF

gen_res_top.py

    A topology is generated for a Tyrosine mimetic for the Amber 14 SB force
    field (ff). A print statement notifies which parameters must be added to 
    the ff (atom types, bonds, angles, (improper) dihedrals, nonbonded).     
    The .pdb file is also edited to label residues as residue, NME or ACE.

CHARGE FITTING PROCEDURE:
1) Obtain alpha & beta constrained simulations

   stage2_sims.sh
		
2) Extract conformations at 1 nm intervals from the simulations

   extract_conformations.py

3) Perform single point HF/6-31G* calculation on each conformation

   single_pt.sh

4) Perform a RESP fitting for each conformation

   resp_fit.sh

5) From the RESP fit charges of all conformations, both alpha and beta, find the average charge
   This step completes stage 2 of the charge parameterization and can be applied to the ff
   
   ave_charges.py
   
    	Finds the average charges from multiconformational RESP fits. 
	The RESP output file from the second iteration, 'resp2.out' is used to retrieve 
	the fitted charge from each configuration.

6) Verify the output with Ramachandran and chi_1 dihedral plots

   ramachandran.py  

