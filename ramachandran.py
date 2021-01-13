#!/usr/bin/env python
# coding: utf-8

import numpy as np
import mdtraj as md
import os
import matplotlib.pyplot as plt

mims = ["FY1a", "FY1b", "FY2a", "FY2b", "FY3a", "FY3b", "FY5a", "FY5b", "FY6a"]# "FY6b"]

home = os.getcwd() 
font = {'color': 'black', 'weight': 'semibold', 'size': 20}

for mim in mims:

    # Enter directory for mimetic
    #os.chdir("/home/finnl92/thesis/stage2chg-param/" + mim)
    os.chdir("C:/Users/Lauren/Documents/Thesis/stage2chg-param/" + mim)
    
    # Obtain the dihedrals from the trajectory
    traj = md.load("fTyr_md.xtc", top="fTyr_md.gro").remove_solvent()
    atoms, bonds = traj.topology.to_dataframe()
    phi_indicies, psi_indicies = [5, 7, 9, 26], [7, 9, 26, 33]
    angles = md.compute_dihedrals(traj, [phi_indicies, psi_indicies])
    
    # Create a plot
    fig, ax = plt.subplots(figsize=(15,10))
    plt.scatter(angles[:, 0], angles[:, 1], marker='x', c=traj.time/1000)
    plt.xlabel(r'$\Phi$ Angle [radians]', fontdict=font, labelpad=10)
    cbar = plt.colorbar()
    cbar.set_label('Time [ns]')
    plt.xlim(-np.pi, np.pi)
    plt.ylabel(r'$\Psi$ Angle [radians]', fontdict=font, labelpad=10)
    plt.ylim(-np.pi, np.pi)
    plt.tick_params(axis='y', labelsize=16, direction='in', width=2, \
                    length=5)
    plt.tick_params(axis='x', labelsize=16, direction='in', width=2, \
                    length=5)
    for i in ["top","bottom","left","right"]:
        ax.spines[i].set_linewidth(3) 
    plt.grid(True)

    plt.savefig("ramachandranplot.png")
    plt.show()
