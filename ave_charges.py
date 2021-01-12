# -*- coding: utf-8 -*-

import os
from sys import argv, exit
from numpy import mean, std
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

def main(argv):
    """
    Finds the average charges from multiconformational RESP fits.

    Subdirectories for the configurations are named 'config_xx'. The .pdb in 
    config_0 directory is used to define the atoms by number in the dictionary
    Residue.atoms. The RESP output file from the second iteration, 'resp2.out' is
    used to retrieve the fitted charge from each configuration.
   
    Parameters
    ----------
    argv[1] : str
        Integer as a string, corresponding to the mimetic's number. E.g. '2' for FY2

    Returns
    -------
    1) Data on multiconformational RESP fit under: 
        ~/thesis/stage2ch-param/stage2_output
    2) Updated charges for ff under ~/thesis/amber14sb.ff/mimetics.rtp
    3) A plot is generated for the charges under:
        ~/thesis/stage2ch-param/stage2_output

    """

    # Check command line arguement
    if len(argv) != 2:
        print("Usage: ave_charges.py <int>")
        exit(1)

    # Which mimetic?
    which = str(argv[1])
    name = "FY" + which

    # Set up path variables
    os.chdir('/home/finnl92/thesis/stage2chg-param')
    alpha = name + "a/configurations"
    beta = name + "b/configurations"
    os.chdir(alpha)
    alpha_home = os.getcwd()
    os.chdir('../../' + beta)
    beta_home = os.getcwd()
    
    # Read in a pdb file
    os.chdir(alpha_home + "/config_0")
    pdb_lines = read_file("geo.pdb")
    for line in pdb_lines:
        line = line.split()
        if "ATOM" in line:
            Residue.pdb.append(line)
    
    # Use the .pdb file to define atoms of the residue
    for line in Residue.pdb[6:-6]:
        Residue.atoms[int(line[1])] = Residue(line[1], line[2])
    os.chdir(alpha_home)
    
    # Collect all charges from alpha conformations for each atom in a list
    folders = [f for f in os.listdir() if "config" in f]
    for folder in folders:
        os.chdir(folder)
        include_charges()
        os.chdir(alpha_home)
        
    # Collect charges from the beta conformations
    os.chdir(beta_home)
    folders = [f for f in os.listdir() if "config" in f]
    for folder in folders:
        os.chdir(folder)
        include_charges()
        os.chdir(beta_home)
    
    # Include any canonical charges    
    canonical_tyr()
    
    # Output info on atoms to file
    data = [str(item) + "\n" for _, item in Residue.atoms.items()]
    os.chdir("/home/finnl92/thesis/stage2chg-param/stage2_output")
    output = open(name, mode="w")
    output.writelines(data)
    output.close()

    # Make plot 
    make_plot(which)
   
    # Edit mimetics.rtp in ff with new charges
    replace_charges(name)

class Residue():
    """
    The atoms belonging to the residue are stored within the dictionary 
    Residue.atoms, where individual atoms are accessed using their integer 
    values from the .pdb file as the key. 
    
    """
    pdb = []
    atoms = dict()
    
    def __init__(self, num, name):
        self.num = int(num)
        self.name = name
        self.charges = []
        self.element = "unknown"
        self.canonical = None

    def __str__(self):
        return "Number: {}, Name: {}, Element: {}, Configs: {}, Charge: {}00, STD: {}"\
            .format(self.num, self.name, self.element, \
                len(self.charges), round(mean(self.charges), 4),\
                    round(std(self.charges), 4))
        
def read_file(file_name):
    """
    Return the lines of a file in a list.
    
    Parameters
    ----------
    file_name : file
        A text file, such as a .pdb or .mol2 file.
    
    Returns
    -------
    result : string list
        An ordered, line-separated list of strings from the file.
    
    """
    try:     
        with open(file_name, mode='r') as file: 
            lines = file.readlines()
    except OSError:
        print("file", file_name, "not found.")
    file_lines = []
    for line in lines:
        file_lines.append(line)
    return file_lines
            
def include_charges():
    """
    Include RESP fit charges for a configuration in Residue.atoms[x].charges
    
    For each atom 'x' in the residue, the RESP fitted charge for the 
    conformation in the current working directory is added to the charges list.
    The second RESP fit should be present in the folder as "resp2.out".
    
    """
    resp_out = read_file("resp2.out")
    
    # find the relevant lines from output for charges
    for i, line in enumerate(resp_out):
        if "Point Charges Before & After Optimization" in line:
            begin = i + 3
        if "Sum over the calculated charges:" in line:
            end = i - 1
    charges = resp_out[begin + 6 : end - 6]
    
    # append the charges to the matching atom's charge list under: 
    # Residue.atoms[x].charges
    check = False
    for line in charges:
        number = int(line.split()[0])
        charge = float(line.split()[3])
        if -1.5 > charge or charge > 1.5: 
            check = True
        Residue.atoms[number].charges.append(charge) 
        Residue.atoms[number].element = line.split()[1]
    if check:
        print("Check RESP output: ", os.getcwd())

def replace_charges(name):
    """
    Replace charges in 'mimetics.rtp' in the ff with the new charges.

    """
    mim_file = "/home/finnl92/thesis/amber14sb.ff/mimetics.rtp"
    mimetics = read_file(mim_file)

    # find the relevant lines
    start = None
    for i, line in enumerate(mimetics):
        if name in line:
            start = i + 2
        if start != None and "bonds" in line:
            end = i
            break
    
    atoms = mimetics[start:end]
    update = []
    for atom in atoms:
        key = int(atom.split()[-1]) + 6
        new_charge = "{:.4f}".format(round(mean(Residue.atoms[key].charges),4))
        new_charge = "{:>{f}}00".format(str(new_charge), f=7)
        update.append(atom[:21] + new_charge + atom[30:])
    mimetics[start:end] = update
    update_mim = open(mim_file, mode="w")
    update_mim.writelines(mimetics)
    update_mim.close()
    return None

def canonical_tyr():
    """
    Include the canonical charge for atoms corresponding to Tyr.
    
    Some of the atoms in the mimetic correspond to an atom on canonical
    Tyrosine. For those atoms, include its canonical charge under the 
    canonical attribute to the Residue class object. 

    """
    
    # Get the charge data for canonical Tyrosine
    rtp_file = "/home/finnl92/thesis/amber14sb.ff/aminoacids.rtp"
    amino_acids = read_file(rtp_file)
    start = None
    for i, line in enumerate(amino_acids):
        if "[ TYR ]" in line:
            start = i + 2
        if start != None and "bonds" in line:
            end = i
            break
    tyrosine = amino_acids[start:end]
    
    # Update the canonical attribute for matching atoms
    for line in tyrosine:
        tyr_name, _, tyr_charge, tyr_num = line.split()
        for i, atom in Residue.atoms.items():
            if atom.name == tyr_name:
                atom.canonical = tyr_charge
                break

def make_plot(num):
    """
    Make a plot of charges, including std error bar and picture of \
    structure.

    Parameters
    ----------
    name : str
        Name of the structure. e.g. "FY1", "FY2" etc.
    
    """
    font = {'color': 'black', 'weight': 'semibold', 'size': 20}
    os.chdir("/home/finnl92/thesis/stage2chg-param/images")
    
    # Data for the mimetic charges
    names = [atom.name for _, atom in Residue.atoms.items()]
    charges = [round(mean(atom.charges), 4) for _, atom in \
            Residue.atoms.items()]

    # Data for canonical tyrosine charges
    stdev = [round(std(atom.charges), 4) for _, atom in Residue.atoms.items()]
    names2 = [atom.name for _, atom in Residue.atoms.items() if \
            atom.canonical != None]
    canonical  = [float(atom.canonical) for _, atom in Residue.atoms.items() \
            if atom.canonical != None]
    
    # Include molecule images
    im1 = "FY" + num + ".png"
    tyr = "tyrosine.png"
    mimetic = mpimg.imread(im1)
    tyrosine = mpimg.imread(tyr)

    # Add figure and subplots
    fig = plt.figure(constrained_layout=True, figsize=(15,10))
    gs = GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])   
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1:, :])

    # Add plots and customize figure
    ax3.scatter(names, charges, linewidth=3, label="FpY " + num, marker='o')
    ax3.scatter(names2, canonical, linewidth=3, label="Tyrosine", color='r', marker='*')
    plt.errorbar(names, charges, yerr=stdev, linewidth=3, fmt="none")
    ax3.tick_params(axis='y', labelsize=16, direction='in', width=2, \
                    length=5)
    ax3.tick_params(axis='x', labelsize=16, direction='in', width=2, \
                    rotation=45, length=5)
    plt.xlabel("Atom Names", fontdict=font, labelpad=10)
    plt.ylabel("Atomic Partial Charge (e)", fontdict=font, labelpad=10)
    plt.legend(prop={"size": 15})
    ax3.grid(True)
    ax1.imshow(mimetic)
    ax2.imshow(tyrosine)
    ax2.text(0,0,"Tyrosine", verticalalignment="bottom", fontdict=font, \
        color="red")
    ax1.text(0,0,"pTyrosine mimetic " + num, verticalalignment="bottom", \
        fontdict=font, color="red")
    for ax in [ax1, ax2]:
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
    for i in ["top","bottom","left","right"]:
        ax3.spines[i].set_linewidth(3) 

    # Save figure and terminate function    
    plt.savefig("FY" + num + "_charges.png")
    return None

if __name__ ==  '__main__':
    main(argv)
    exit(0)
