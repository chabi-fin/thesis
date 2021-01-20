#!/usr/bin/env python
# coding: utf-8

import numpy as np
import mdtraj as md
import os
import matplotlib.pyplot as plt

# Dictionary for phi, psi and chi_1 indicies
mims = {"FY1" : ((5,7,9,26),(7,9,26,33),(7,9,11,14)), "FY2" : ((5,7,9,26), \
        (7,9,26,33),(7,9,11,14)), "FY3" : ((5,7,9,26),(7,9,26,35), \
        (7,9,11,14)), "FY5" : ((5,7,9,27),(7,9,27,41),(7,9,11,14)), \
        "FY6" : ((5,7,9,26),(7,9,26,32),(7,9,11,14))}

# Set up variables
home = os.getcwd() 
font = {'color': 'black', 'weight': 'semibold', 'size': 20}
conformations = []
for mim in mims.keys():
    conformations.append(mim + "a")
    conformations.append(mim + "b")

# Dihedral angles for both conformations of each mimetic
for mim in conformations:
    
    # Mimetic name, e.g. "FY1"
    base = mim[:-1]
    
    # Enter directory for mimetic, top : qcm, bottom : local
    os.chdir("/home/finnl92/thesis/stage2chg-param/" + mim)
    #os.chdir("C:/Users/Lauren/Documents/Thesis/stage2chg-param/" + mim)
    
    # Obtain the dihedrals from the trajectory
    traj = md.load("fTyr_md.xtc", top="fTyr_md.gro").remove_solvent()
    atoms, bonds = traj.topology.to_dataframe()
    phi_indicies, psi_indicies, chi_indicies = mims[base]
    angles = md.compute_dihedrals(traj, [phi_indicies, psi_indicies, \
                                         chi_indicies])
    
    # Get the position restraint vector
    with open("posre.itp", mode="r") as file:
        lines = file.readlines()
    for i, line in enumerate(lines):
        if i == 4: 
            vector = line.split()[2:]
    force_const = " ".join(vector) + r" kJ/mol$\cdot$nm$^{2}$"

    # Create a Ramachandran plot
    fig, ax = plt.subplots(figsize=(15,10))
    plt.scatter(angles[:, 0], angles[:, 1], marker='x', c=traj.time/1000)
    plt.xlabel(r'$\Phi$ Angle [radians]', fontdict=font, labelpad=10)
    cbar = plt.colorbar()
    cbar.set_label('Time [ns]', fontdict=font, labelpad=10)
    cbar.ax.set_yticklabels(np.arange(0,100,20),fontsize=16)
    cbar.ax.tick_params(labelsize=16, direction='out', width=2, length=5)
    cbar.outline.set_linewidth(3)
    plt.xlim(-np.pi, np.pi)
    plt.ylabel(r'$\Psi$ Angle [radians]', fontdict=font, labelpad=10)
    plt.ylim(-np.pi, np.pi)
    plt.tick_params(axis='y', labelsize=20, direction='in', width=2, \
                    length=5, pad=10)
    plt.tick_params(axis='x', labelsize=20, direction='in', width=2, \
                    length=5, pad=10)
    for i in ["top","bottom","left","right"]:
        ax.spines[i].set_linewidth(3) 
    plt.grid(True)
    plt.text(0.95, 0.05, "Force constant vector\n" + force_const, \
            fontdict={"size" : 16}, transform=ax.transAxes, \
            bbox=dict(boxstyle="round", fc="gainsboro", alpha=0.5), ma ="center", va="bottom", ha="right")
    if mim[3] == "a":
        conf = "alpha"
    else:
        conf = "beta"
    plt.title("Fluorinated pTyr " + mim[2] + ", Confomation: " + conf, fontdict=font)
    plt.savefig("ramachandranplot_" + mim[2:] + ".png")
    plt.show()
    
    # Create Chi_1 plot
    fig, ax = plt.subplots(figsize=(15,10))
    plt.hist(angles[:, 2], bins=75, color="darkgoldenrod")
    plt.xlabel(r'$\chi_1$ Angle [radians]', fontdict=font, labelpad=10)
    plt.xlim(-np.pi, np.pi)
    plt.title("Fluorinated pTyr " + mim[2] + ", Confomation: " + conf, fontdict=font)
    plt.tick_params(axis='y', labelsize=20, direction='in', width=2, \
                    length=5, pad=10)
    plt.tick_params(axis='x', labelsize=20, direction='in', width=2, \
                    length=5, pad=10)
    for i in ["top","bottom","left","right"]:
        ax.spines[i].set_linewidth(3) 
    plt.grid(True)
    plt.savefig("chi1_plot_" + mim[2:] + ".png")
