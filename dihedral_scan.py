# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 12:20:53 2021

@author: Lauren
"""

import mytools
import matplotlib.pyplot as plt
import numpy as np
from sys import argv, exit

# =============================================================================
# #Assign varaible names to command line arguements
# if len(sys.argv) != 2:
#     print("Usage: python extract_scan_summary.py $PATH_FOR_SCAN.LOG")
#     sys.exit(1)
# path = sys.argv[1]
# =============================================================================

path = "type3.log"

# load log file into the memory
log_file = mytools.read_file(path)

# Grab input settings, scan_incs : # of scan increments,
# spacing : degrees seperation for each step in scan,
# initial_value : intial value in (deg) for the scan 
for i, line in enumerate(log_file):
    if "The following ModRedundant input section has been read:" in line:
        input_settings = log_file[i + 1].split()
    if "!    Initial Parameters    !" in line:
        initial_param = i
        break
for line in log_file[i:]:
    if "Scan" in line:
        initial_value = line.split("Scan")[0].split()[-1]
        break
scan_incs = int(input_settings[6])
spacing = float(input_settings[7])
initial_value = float(initial_value)

# Extract the Eigenvalues into : list(eigenvalues). 
# base_energy stores the energy which should be added to all energies
for i, line in enumerate(log_file):
    if "Summary of Optimized Potential Surface Scan" in line:
        base_energy = line.split("add")[1]
        base_energy = float(base_energy.split()[0])
        start_summary = i
        break
eigenvalues, j = [], 0
while ("GradGradGradGrad" not in log_file[start_summary + j]):
    line = log_file[start_summary + j]
    if "Eigenvalues" in line:
        eigens = line.split()[2:]
        eigenvalues.extend(eigens)
    j += 1
eigenvalues = list(map(float, eigenvalues))

# Create list(degrees) the degree value at each scan point
degrees = [initial_value]
current_val = initial_value
for i in range(scan_incs):
    current_val += spacing
    degrees.append(current_val)

# Generate a plot for the scan
# Convert Hartree to kcal/mol
degrees = np.array(degrees[1:]) - 360
eigenvalues = (np.array(eigenvalues[1:]) - min(eigenvalues)) * 627.503
plt.plot(degrees, eigenvalues)
plt.grid()
plt.xlabel("degrees")
plt.ylabel("potential (kcal/mol)")

# Print the potential as difference b/w min and max energy
amplitude =(max(eigenvalues)-min(eigenvalues))
print("V_n : {}".format(amplitude))
