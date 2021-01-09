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
    
    Parameters
    ----------
    Subdirectories for the configurations are named 'config_xx'. The .pdb in 
    config_0 directory is used to define the atoms by number in the dictionary
    Residue.atoms. The RESP output file from the second iteration, 'resp2.out' is
    used to retrieve the fitted charge from each configuration.
    
    Returns
    -------
    result : 
    
    """
    
    # Which mimetic?
    which = "2"
    
    # Set up path variables
    alpha = "FY" + which + "a/configurations"
    beta = "FY" + which + "b/configurations"
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
    os.chdir("../..")
    
    # Print out info on atom, including its average charge
    print(*[str(item) for _, item in Residue.atoms.items()], sep="\n")    
    
    # Make plot 
    make_plot("FY" + which)
    
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

def make_plot(name):
    """
    Make a plot of charges, including std error bar and picture of \
    structure.

    Parameters
    ----------
    name : str
        Name of the structure. e.g. "FY1", "FY2" etc.
    
    """
    font = {'family': 'Times New Roman', 'color': 'black', 'weight': \
            'semibold', 'size': 20}
        
    os.chdir("images")

    names = [atom.name for _, atom in Residue.atoms.items()]
    charges = [round(mean(atom.charges), 4) for _, atom in Residue.atoms.items()]
    stdev = [round(std(atom.charges), 4) for _, atom in Residue.atoms.items()]
    im = name + ".png"
    molecule = mpimg.imread(im)
    
    fig = plt.figure(constrained_layout=True, figsize=(15,10))
    
    gs = GridSpec(3, 1, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])   

    #plt.figure(figsize=(10,7))
    ax2 = fig.add_subplot(gs[1:, 0])
    #ax2 = plt.axes()
    ax2.scatter(names, charges, label="FpY" + name[2])
    plt.errorbar(names, charges, yerr=stdev, fmt="none")
    ax2.tick_params(axis='y', labelsize=16, direction='in', width=2, \
                    length=5)
    ax2.tick_params(axis='x', labelsize=16, direction='in', width=2, \
                    rotation=45, length=5)
    plt.xlabel("Atom Names", fontdict=font, labelpad=10)
    plt.ylabel("Atomic Partial Charge (e)", fontdict=font, labelpad=10)
    plt.legend()
    ax1.imshow(molecule)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    
    for i in ["top","bottom","left","right"]:
        ax2.spines[i].set_linewidth(2) 
    plt.savefig(name + "_charges.png")
    #plt.show()

if __name__ ==  '__main__':
    main(argv)
    exit(0)
