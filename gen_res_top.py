# -*- coding: utf-8 -*-
""" 
Created on Tue Oct  6 20:23:49 2020

@author: Lauren
"""    

from sys import argv, exit
import string

def main(argv):
    """
    Generate a residue topology and print missing parameters.
    
    A topology is generated for a Tyrosin mimetic for the Amber 14 SB force
    field (ff). A print statement notifies which parameters must be added to 
    the ff (atom types, bonds, angles, (improper) dihedrals, nonbonded).     
    The .pdb file is also edited to label residues as residue, NME or ACE.
    
    Parameters
    ----------
    file_itp : file
        A .itp file, e.g. name_GMX.itp.
    
    residue : str
        3 or 4 letter string to label the residue.
        
    bonded_itp : file
        A .itp file for bonded ff parameters.
    
    nonbonded_itp : file
        A .itp file for nonbonded ff parameters.
    
    Returns
    -------
    result : print statements, edits .pdf file
        Prints the residue topology for the mimetic.
        Prints the parameters which must be added to the forcefield.
        Edits the .pdb file to update residue names. 
    
    """
    # if len(argv) != 3:
    #     print('Usage: residue.py <itp file> <residue name>')
    #     exit(1)
    # basename = argv[1]
    # ffbonded = argv[2]
    # ffnonbonded = argv[3]
    # global residue_name
    # residue_name = argv[4].upper()
    # if len(residue_name) not in [3, 4]:
    #     print("Residue name should have a length of 3 or 4.")
    #     exit(2)
    
    ##### Values for testing #####
    basename = "FY5"
    folder = basename + "/" + basename + ".acpype/"
    file_itp = folder + basename + "_GMX.itp" 
    file_pdb = folder + basename + "_NEW.pdb"
    ffbonded = "../amber14sb.ff/ffbonded.itp"
    ffnonbonded = "../amber14sb.ff/ffnonbonded.itp"
    global residue_name
    residue_name = "FY5"
    ##############################
    
    # Read the .itp file into Atom.itp
    lines = read_file(file_itp)
    for line in lines:
        Atom.itp.append(line)
    blocks_dict() # store block indicies in itp file : "Atom.blocks"    
        
    # initialize all atoms in the Atom class
    start, end = Atom.blocks["atoms"]
    for atom in Atom.itp[start : end + 1]:
        Atom.atoms.append(Atom(atom.split()[0]))
    for atom in Atom.atoms:
        Atom.original_names[atom.original] = atom
    
    # Change the residue names for the caps
    res_to_NME()
    res_to_ACE()
    
    # Assign Tyrosine-like atom names in the residue
    Tyrosine.assign_tyr_types()
    mimetic_atomnames()
    
    # Print the parameters that need to be added to the force field
    atomtypes_dict()
    check_bonds(ffbonded)
    check_angles(ffbonded)
    check_dihedrals(ffbonded)
    
    for atom in Atom.atoms:
        print(atom)
    
    reorder_atoms()
    update_pdb(file_pdb)
    
    # Print the residue topology which can be added to aminoacids.rtp
    residue_rtp = ["[ {} ]\n".format(residue_name)] + Mimetic.atoms() + \
        Mimetic.bonds() + Mimetic.impropers()
    print(*residue_rtp, sep="")
    
    return None
    
class Atom:
    """
    An instance in the Atom class, represents an atom in the .itp file.
    
    Attributes
    ----------
    self.num : int
        The number assigned to the atom in the .itp file.
    self.element : str
        The element symbol.
    self.itp_index : int
        The line index for the atom in the .itp file, under [ atoms ].

    """
    
    itp = []
    gro = []
    atoms = []
    original_names = dict()
    blocks = dict()
    gaff_amber_types = {"ca" : "CA", "c3" : "CT", "p5" : "P", "f" : "F", \
                        "hc" : "HC", "h1" : "H1", "o" : "O", "ha" : "HA" }
    residue_names = set()
    reordered = []
    
    def __init__(self, num):
        self.num = int(num)
        self.itp_index = Atom.blocks["atoms"][0] + self.num - 1
        self.itp_line = Atom.itp[self.itp_index].split()
        # element is the atom with the numeric identifier removed
        self.element = self.itp_line[4].rstrip(string.digits)
        self.atomtype = self.itp_line[1]
        self.name = self.itp_line[4].lower()
        self.charge = self.itp_line[6]
        self.residue = residue_name
        self.original = self.name.upper()
        
    def __str__(self):
        return "Number: {}, Element: {}, Name: {}, Type: {}, Residue: {}".\
            format(self.num, self.element, self.name, self.atomtype, self.residue)
        
    def bonds(self):
        '''list of atom numbers the atom is bound to.'''
        i = Atom.blocks["bonds"][0] # start index for "[ bonds ]"
        j = Atom.blocks["bonds"][1] # end index for "[ bonds ]"
        links = []
        for bond in Atom.itp[i : j + 1]:
            pair = bond.split()[:2]
            if str(self.num) in pair:
                for a in pair:
                    if a != str(self.num):
                        links.append(a)
        return links

    def neighbors(self):
        """list of atoms the atom is bound to."""
        nums = list(map(lambda x : Atom.atoms[int(x) - 1], self.bonds()))
        return nums
    
    def is_methyl(self):
        '''true if atom is the carbon of a methyl group'''
        if self.element != 'C':
            return False
        count_h = 0
        for n in self.neighbors():
            if n.element == 'H':
                count_h += 1
        if count_h == 3:
            return True
        return False
    
    def is_carbonyl(self):
        '''true if atom is the carbon of a carbonyl group'''
        if self.element != 'C':
            return False
        ns = self.neighbors()
        if len(ns) != 3:
            return False
        for n in ns:
            if n.element == "O" and len(n.bonds()) == 1:
                return True
        return False
     
class Mimetic:
    """
    """
    
    def atoms():
        """
        Generate the " [ atoms ]" section for the residue.
        
        Returns
        -------
        result : str list
            Returns the atoms section for the residue, to be added to the 
            aminoacids.rtp file.
            
        """
        atoms = [" [ atoms ]\n"]
        for i, atom in enumerate(Atom.reordered):
            atom_name = "{:>{field_size}}".format(atom.name, field_size = 6)
            atom_type = "{:<{field_size}}".format(atom.atomtype, field_size = 2)
            charge = "{:>{field_size}}".format(atom.charge, field_size = 18)
            atom_num = "{:>{field_size}}".format(str(i + 1), field_size = 6)
            line = atom_name + "    " + atom_type + charge + atom_num + "\n"
            atoms.append(line)
        return atoms

    def bonds():
        """
        Generate the " [ bonds ]" section for the residue.
        
        Returns
        -------
        result : str list
            Returns the bonds section for the residue, to be added to the 
            aminoacids.rtp file.

        """
        ### Need bonds connecting to backbone!!!!!!
        bonds = [" [ bonds ]\n"]
        pairs = []
        for atom in Atom.reordered: 
            for n in atom.neighbors():
                if n.residue == residue_name:
                    pair = sorted((atom.name, n.name))
                    if pair not in pairs:
                        pairs.append(pair)
        for pair in pairs:
            first = "{:>{field_size}}".format(pair[0], field_size=6)
            second = "{:>{field_size}}".format(pair[1], field_size=6) 
            bonds.append(first + second + "\n")
        bonds.append("    -C     N\n")
        return bonds

    def impropers():
        """
        Generate the " [ impropers ]" section for the residue.

        result : str list
            Returns the improper dihedrals section for the residue, to be added to the 
            aminoacids.rtp file.

        """
        impropers = [" [ impropers ]\n"] #, "    -C    CA     N     H\n", \
                     #"    CA    +N     C     O\n"]
        start, end = Atom.blocks["improper dihedrals"]
        for line in Atom.itp[start : end + 1]:
            improper = line.split(";")[1].split()
            improper = [Atom.original_names[x.strip("-")] for x in improper]
            improper = [i.name for i in improper]
            improper = "{:>{f}}{:>{f}}{:>{f}}{:>{f}}\n".format(*improper, f=6)
            if "C    CA     N     H" in improper:
                improper = "    -C    CA     N     H\n"
            elif "CA     N     C     O" in improper:
                improper = "    CA    +N     C     O\n"
            impropers.append(improper)
        return impropers    
    
class Tyrosine:
    
    def tyr_CA():
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                return None
            if atom.atomtype == "c3":
                conds_cx = 4 * [False]
                for n in atom.neighbors():
                    if n.is_carbonyl():
                        conds_cx[0] = True
                    if n.atomtype == "h1" or "H1":
                        conds_cx[1] = True
                    if n.atomtype == "n" or "N":
                        conds_cx[2] = True
                    if n.atomtype == "c3" or "CT" or "CX": 
                        conds_cx[3] = True
                if all(conds_cx):
                    atom.atomtype = "CX"
                    atom.name = "CA"
                    return True
        return False
                
    def tyr_HA():
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.atomtype == "h1":
                conds_ha = 2 * [False]
                if len(atom.bonds()) == 1:
                    conds_ha[0] = True
                if atom.neighbors()[0].name == "CA":
                    conds_ha[1] = True
                if all(conds_ha):
                    atom.atomtype = "H1"
                    atom.name = "HA"
                    return True
        return False
                    
    def tyr_CB():
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.atomtype == "c3":
                conds_ct = 4 * [False]
                for n in atom.neighbors():
                    if n.name == "CA":
                        conds_ct[0] = True
                    if n.atomtype == "ca":
                        conds_ct[1] = True
                    if n.atomtype == "hc":
                        conds_ct[2] = True
                    if n.atomtype == "hc" and conds_ct[2]:
                        conds_ct[3] = True
                if all(conds_ct):
                    atom.atomtype = "CT"
                    atom.name = "CB"
                    return True
        return False
    
    def tyr_HBs():
        first = True
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "H":
                conds_hb = 2 * [False]
                links = atom.neighbors()
                if len(links) == 1:
                    conds_hb[0] = True
                if links[0].name == "CB":
                    conds_hb[1] = True
                if all(conds_hb) and first:
                    atom.atomtype = "HC"
                    atom.name = "HB1"
                    first = False
                elif all(conds_hb):
                    atom.atomtype = "HC"
                    atom.name = "HB2"
                    return True
        return False
    
    def tyr_C(): 
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.is_carbonyl():
                for n in atom.neighbors():
                    if n.name == "CA":
                        atom.atomtype = "C"
                        atom.name = "C"
                        return True
        return False
                        
    def tyr_O():
        for atom in Atom.atoms:
            if atom.atomtype.isupper(): 
                continue
            if atom.element == "O" and atom.neighbors()[0].name == "C":
                atom.atomtype = "O"
                atom.name = "O"
                return True
        return False

    def tyr_N():
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "N":
                for n in atom.neighbors():
                    if n.name == "CA":
                        atom.atomtype = "N"
                        atom.name = "N"
                        return True
        return False
    
    def tyr_H():
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "H" and atom.neighbors()[0].name == "N":
                atom.atomtype = "H"
                atom.name = "H"
                return True
        return False
    
    def tyr_CG():
        for atom in Atom.atoms:
            if atom.atomtype.isupper(): 
                continue
            if atom.atomtype == "ca":
                conds_cg = 3 * [False]
                first = True
                for n in atom.neighbors():
                    if n.name == "CB":
                        conds_cg[0] = True
                    if n.atomtype == "ca" and first:
                        conds_cg[1] = True
                        first = False
                    elif n.atomtype == "ca" and (not first):
                        conds_cg[2] = True
                if all(conds_cg):
                    atom.atomtype = "CA"
                    atom.name = "CG"
                    return True
        return False    
    
    def tyr_CDs():
        first = True
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "C":
                conds_cd = 3 * [False]
                for n in atom.neighbors():
                    if n.element == "H":
                        conds_cd[0] = True
                    if n.name == "CG":
                        conds_cd[1] = True
                    if n.atomtype == "ca":
                        conds_cd[2] = True
                if all(conds_cd) and first:
                    atom.atomtype = "CA"
                    atom.name = "CD1"
                    first = False
                elif all(conds_cd):
                    atom.atomtype = "CA"
                    atom.name = "CD2"
                    return True
        return False
    
    def tyr_CEs():
        both = 0
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "C":
                conds_ce = 2 * [False]
                for n in atom.neighbors():
                    if n.atomtype == "ca":
                        conds_ce[0] = True
                    elif n.element == "H":
                        conds_ce[1] = True
                    else: 
                        name = n.name
                    if all(conds_ce) and name == "CD1":
                        atom.atomtype = "CA"
                        atom.name = "CE1"
                        both += 1
                    elif all(conds_ce) and name == "CD2":
                        atom.atomtype = "CA"
                        atom.name = "CE2"
                        both += 1
        if both == 2:
            return True
        return False
    
    def tyr_HCs():
        all = 0
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "H":
                n = atom.neighbors()
                name = n[0].name
                cas = ["CD1", "CD2", "CE1", "CE2"]
                has = ["HD1", "HD2", "HE1", "HE2"]
                if len(n) == 1:
                    cond = True
                if name in cas and cond:
                    atom.atomtype = "HA"
                    for c, h in zip(cas, has):
                        if c == name:
                            atom.name = h
                    all += 1
        if all == 4:
            return True
        return False
    
    def tyr_CZ():
        second = False
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "C":
                conds_cz = 4 * [False]
                for n in atom.neighbors():
                    if n.atomtype == "c":
                        conds_cz[0] = True
                    elif n.name is "CE1" or "CE2":
                        if second == True:
                            conds_cz[1] = True
                        conds_cz[2] = True
                        second = True
                    elif n.element == "O":
                        conds_cz[3] = True
                if all(conds_cz):
                    atom.atomtype = "C"
                    atom.name = "CZ"
                    return True
        return False

    def tyr_OH():
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "O":
                conds_oh = 2 * [False]
                for n in atom.neighbors(): 
                    if n.name == "CZ":
                        conds_oh[0] = True
                    elif n.element == "H":
                        conds_oh = True
                    if all(conds_oh):
                        atom.name = "OH"
                        atom.atomtype = "OH"
                        return True
        return False
    
    def tyr_HH():
        for atom in Atom.atoms:
            if atom.atomtype.isupper():
                continue
            if atom.element == "H" and atom.neighbors()[0].name == "OH":
                atom.name = "HH"
                atom.atomtype = "HO"
                return True
        return False
                
    def assign_tyr_types():
        Tyrosine.tyr_CA()
        Tyrosine.tyr_HA()
        Tyrosine.tyr_CB()
        Tyrosine.tyr_HBs()
        Tyrosine.tyr_C()
        Tyrosine.tyr_O()
        Tyrosine.tyr_N()
        Tyrosine.tyr_H()
        Tyrosine.tyr_CG()
        Tyrosine.tyr_CDs()
        Tyrosine.tyr_CEs()
        Tyrosine.tyr_HCs()
        Tyrosine.tyr_CZ()
        Tyrosine.tyr_OH()
        Tyrosine.tyr_HH()

def read_file(file_name):
    """
    Return the lines of a (.pdb or .mol2) file in a list.
    
    Parameters
    ----------
    file_name : file
        A .pdb or .mol2 file.
    
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

def blocks_dict():
    """
    Create a dictionry for indicies in .itp file. Stored under "Atom.itp".
    
    Returns
    -------
    result : dictionary with list values
        Structure: Atom.blocks["[ key ]"] = [<index start data>, <index end data>]
        
    """
    index = 0 # keep track of index
    for line in Atom.itp:
        if "[ atomtypes ]" in line:
            end = index
            while Atom.itp[end].split() != []:
                end += 1
            Atom.blocks["atomtypes"] = [index + 2, end - 1]
        elif "[ moleculetype ]" in line:
            end = index
            while Atom.itp[end].split() != []:
                end += 1
            Atom.blocks["moleculetype"] = [index + 2, end - 1]
        elif "[ atoms ]" in line:
            end = index
            while Atom.itp[end].split() != []:
                end += 1
            Atom.blocks["atoms"] = [index + 2, end - 1]
        elif "[ bonds ]" in line:
            end = index
            while Atom.itp[end].split() != []:
                end += 1
            Atom.blocks["bonds"] = [index + 2, end - 1]
        elif "[ pairs ]" in line:
            end = index
            while Atom.itp[end].split() != []:
                end += 1
            Atom.blocks["pairs"] = [index + 2, end - 1]
        elif "[ angles ]" in line:
            end = index
            while Atom.itp[end].split() != []:
                end += 1
            Atom.blocks["angles"] = [index + 2, end - 1]
        elif "[ dihedrals ] ; propers" in line:
            end = index
            while Atom.itp[end].split() != []:
                end += 1
            Atom.blocks["proper dihedrals"] = [index + 3, end - 1]
        elif "[ dihedrals ] ; impropers" in line:
            end = len(Atom.itp) - 1
            Atom.blocks["improper dihedrals"] = [index + 3, end]            
        index += 1

def find_NME():
    """
    Return the nitrogen atom from the NME cap.
    
    Returns
    -------
    result : Atom
        The N atom from the NME residue. See class Atom.
    
    """
    nitrogens = []
    for atom in Atom.atoms:
        if atom.element == 'N':
            nitrogens.append(atom)
    for n in nitrogens:
        conds = 3 * [False]
        for bond in n.bonds():
            if Atom(bond).is_methyl():
                conds[0] = True
            if Atom(bond).element == 'H':
                conds[1] = True
            if Atom(bond).is_carbonyl():
                conds[2] = True
        if all(conds):
            return n
    return None
        
def find_ACE():
    """
    Return the carbon atom from the ACE cap.
    
    Returns
    -------
    result : Atom
        The N atom from the NME residue. See class Atom.
    
    """
    carbonyls = []
    for atom in Atom.atoms:
        if atom.is_carbonyl():
            carbonyls.append(atom)
    for c in carbonyls:
        conds = 2 * [False]
        for bond in c.bonds():
            if Atom(bond).element == 'N':
                conds[0] = True
            if Atom(bond).is_methyl():
                conds[1] = True
        if all(conds):
            return c    
    return None

def res_to_NME():
    """
    Change residue name to NME for atoms in the NME cap.
    
    Returns
    -------
    result : None
        Changes the residue name to "NME" for atoms (under "[ atoms ]" in the 
        .itp file) belonging to the N-methyl amide cap.
    
    """
    nme = find_NME()
    if nme == None:
        print("Error: Could not find NME residue.")
        return False
    nme.label_ind = 0
    changes = [nme]
    names = ("N", "H", "CH3", "HH31", "HH32", "HH33")
    types = ("N", "H", "CT", "H1", "H1", "H1")
    for n in nme.neighbors():
        if n.is_carbonyl():
            continue
        if n.element == 'H':
            changes.append(n)
            n.label_ind = 1
        if n.is_methyl():
            changes.append(n)
            n.label_ind = 2
            x = 0
            for h in n.neighbors():
                if h.element == 'H':
                    changes.append(h)
                    h.label_ind = 3 + x
                    x += 1
    for change in changes:
        i = change.label_ind
        change.residue = "NME"
        change.name = names[i]
        change.atomtype = types[i]
    return True
        
def res_to_ACE():
    """
    Change residue name to NME for atoms in the NME cap.
    
    Returns
    -------
    result : None
        Changes the residue name to "ACE" for atoms (under "[ atoms ]" in the 
        .itp file) belonging to the methyl-acetyl cap.
    
    """
    ace = find_ACE()
    if ace == None:
        print("Error: Could not find ACE residue.")
        return False
    ace.label_ind = 0
    changes = [ace]
    names = ("C", "O", "CH3", "HH31", "HH32", "HH33")
    types = ("C", "O", "CT", "HC", "HC", "HC")
    for n in ace.neighbors():
        if n.element == 'N':
            continue
        if n.element == 'O':
            changes.append(n)
            n.label_ind = 1
        if n.is_methyl():
            changes.append(n)
            n.label_ind = 2
            x = 0
            for h in n.neighbors():
                if h.element == 'H':
                    changes.append(h)
                    h.label_ind = 3 + x
                    x += 1
    for change in changes:
        i = change.label_ind
        change.residue = "ACE"
        change.name = names[i]
        change.atomtype = types[i]

def mimetic_atomnames():
    """
    Assign atomtypes and names to the remaining atoms on the residue.
    
    Returns
    -------
    None.

    """
    total_mims = 0
    for atom in Atom.atoms:
        if atom.residue == residue_name:
            total_mims += 1
        if atom.name.isupper() and atom.residue == residue_name:
            Atom.residue_names.add(atom.name)
    for atom in Atom.atoms:
        if atom.name.isupper():
            continue
        atom.atomtype = Atom.gaff_amber_types[atom.atomtype]
        new_name = atom.element 
        count = 1
        while (new_name in Atom.residue_names):
            strip = "".join([x for x in new_name if x.isalpha()])
            new_name = strip + str(count)
            count += 1
        atom.name = new_name
        Atom.residue_names.add(atom.name)
    if total_mims != len(Atom.residue_names):
        print("\nERROR: NOT ALL ATOMS HAVE UNIQUE NAMES!\n")
    
def atomtypes_dict():
    """
    Check whether atom types should be added to the atomtype dictionary.
    
    Returns
    -------
    result : None
        Checks the atomtypes in the .itp file. Prints a warning and exits the 
        program if there is not already a corresponding Amber atomtype listed 
        in the dictionary Atom.atomtypes.
    
    """
    for atom in Atom.atoms:
        if atom.atomtype.islower():
            print("Please add \"{}\" to the atom type dictionary, "\
                  "Atom.atom_types\n".format(atom.atomtype))
            exit(3)
 
def update_pdb(file):
    """
    Update the .pdb file to include changes to the residues and atom names.

    Parameters
    ----------
    file : .pdb file
        The .pdb file for the residue.

    Returns
    -------
    The .pdb file is updated to include the relevant atom names and residues
    for running a simulation on the residue in the edited amber 14 sb ff.

    """
    
    lines = read_file(file)
    new_pdb = file.strip("_NEW.pdb") + ".pdb"
    file_new = open(new_pdb, mode="w")
    file_new.writelines(lines[0])
    ordered = []
        
    for line in lines[1:-1]:
        old_name = line.split()[2]
        index = int(line.split()[1]) - 1
        old_res = line.split()[3]
        new_res = "{:<{field_size}}".format(Atom.atoms[index].residue, \
                                            field_size=4)
        new_name = "{:<{field_size}}".format(Atom.atoms[index].name, \
                                            field_size=4)
        new_line = line.replace(old_res + (4 - len(old_res))*" ", new_res)
        new_line = new_line.replace(old_name + (4 - len(old_name))*" ", new_name)
        if " " not in new_name:
            new_line = new_line.replace(" " + new_name, new_name + " ")
        ordered.append(new_line)
        #file_new.writelines(new_line)
    ordered.sort(key=lambda x : x.split()[3])
    reorder = ordered
    
    for line in ordered[6:-6]:
        name = line.split()[2]
        for atom in Atom.reordered:
            if atom.name == name:
                index = atom.num
        reorder[index] = line
        
    for i, line in enumerate(reorder):
        line = "ATOM     " + "{:>2}".format(str(i + 1)) + line[11:]
        file_new.writelines(line)

    #file_new.writelines(lines[-1])
    file_new.close()             

def check_bonds(ffbonded):
    """
    Check if all bonds are in the ff and print any missing bonds. 

    Parameters
    ----------
    ffbonded : .itp file
        A copy of the .itp for the bonded parameters in the ff.

    Returns
    -------
    None
        Prints out the necessary edits to the bonds section of the bonded 
        part of the ff.
        
    """
    start, end = Atom.blocks["bonds"]
    pairs = set()
    for line in Atom.itp[start : end + 1]:
        pair = line.split(";")[1].split()
        pair_params = tuple(line.split(";")[0].split()[-2:])
        pair.remove("-")
        type1 = Atom.original_names[pair[0]].atomtype
        type2 = Atom.original_names[pair[1]].atomtype
        pairs.add((tuple(sorted([type1, type2])), pair_params))
        
    lines = read_file(ffbonded)
    bondtypes = []
    for i, line in enumerate(lines):
        if "[ bondtypes ]" in line:
            while lines[i].split() != []:
                bondtypes.append(tuple(sorted(lines[i].split()[:2])))
                i += 1
            break
        
    missing_pairs = []
    for pair in pairs: 
        if pair[0] not in bondtypes:
            b0 = pair[1][0].split("e")
            b0 = round(float(b0[0]) * 10 ** int(b0[1]), 5)
            kb = pair[1][1].split("e")
            kb = round(float(kb[0]) * 10 ** int(kb[1]), 1)
            missing = "  {:<3}{:<11}1    {:<07}   {}".format(pair[0][0], pair[0][1], b0, kb)
            missing_pairs.append(missing)
            
    comment = "  ; added by Finn\n"
    if missing_pairs != []:
        print("The pair(s): ")
        print("  i  j    func       b0        kb")
        print(*missing_pairs, sep=comment, end=comment)
        print("are missing from the forcefield.")
        print("Add them to ffbonded under bondtypes.\n")
    
def check_angles(ffbonded):
    """
    Check if all bond angles are in the ff and print any missing angles. 

    Parameters
    ----------
    ffbonded : .itp file
        A copy of the .itp for the bonded parameters in the ff.

    Returns
    -------
    None
        Prints out the necessary edits to the angles section of the bonded 
        part of the ff.
    
    """
    start, end = Atom.blocks["angles"]
    angles = set()
    for line in Atom.itp[start : end + 1]:
        angle = line.split(";")[1].split()
        angle_params = tuple(line.split(";")[0].split()[-2:])
        angle = [x for x in angle if x != "-"]
        angle = [Atom.original_names[angle[i]].atomtype for i in range(3)]
        angles.add((tuple(angle), angle_params))
    angles = list(angles)
    
    lines = read_file(ffbonded)
    angletypes = []
    for i, line in enumerate(lines):
        if "[ angletypes ]" in line:
            while lines[i].split() != []:
                tup = tuple(lines[i].split()[:3])
                angletypes.append(tup)
                i += 1
            break
        
    missing_angles = []
    for angle in angles:
        if angle[0] not in angletypes and angle[0][::-1] not in angletypes:
            th0 = angle[1][0].split("e")
            th0 = round(float(th0[0]) * 10 ** int(th0[1]), 3)
            cth = angle[1][1].split("e")
            cth = round(float(cth[0]) * 10 ** int(cth[1]), 3)
            ang_atoms = "{:<4}{:<4}{:<13}".format(*angle[0])
            missing = ang_atoms + "1   {:<07}     {:<07}".format(th0, cth)
            missing_angles.append(missing)
              
    comment = "  ; added by Finn\n"
    if missing_angles != []:
        print("The angle(s): ")
        print("i   j   k    func         th0         cth")
        print(*missing_angles, sep=comment, end=comment)
        print("are missing from the forcefield.")
        print("Add them to ffbonded under angletypes.\n")
        
def check_dihedrals(ffbonded):
    """
    Check if all dihedrals are in the ff and print any missing dihedrals. 

    Parameters
    ----------
    ffbonded : .itp file
        A copy of the .itp for the bonded parameters in the ff.

    Returns
    -------
    Prints dihedral tuples which must be added into the .itp file.

    """
    dihedrals = set()
    start, end = Atom.blocks["proper dihedrals"]
    for line in Atom.itp[start : end + 1]:
        dihedral = line.split(";")[1].split()
        hedral_params = tuple(line.split(";")[0].split()[-4:])
        dihedral = [x.strip("-") for x in dihedral]
        dihedral = [Atom.original_names[x].atomtype for x in dihedral]
        dihedrals.add((tuple(dihedral), hedral_params))
    dihedrals = list(dihedrals)
    
    dihedrals_improper = set()
    start, end = Atom.blocks["improper dihedrals"]
    for line in Atom.itp[start : end + 1]:
        dihedral_im = line.split(";")[1].split()
        improper_params = tuple(line.split(";")[0].split()[-4:])
        dihedral_im = [Atom.original_names[x.strip("-")] for x in dihedral_im]
        dihedral_im = [i.atomtype for i in dihedral_im]
        dihedrals_improper.add((tuple(dihedral_im), improper_params))
    
    lines = read_file(ffbonded)
    dihedraltypes = []
    for i, line in enumerate(lines):
        if "[ dihedraltypes ]" in line and "improper" not in line:
            while i < len(lines):
                dihedraltypes.append(tuple(lines[i].split()[:4]))
                i += 1
            break
        
    impropertypes = []
    for i, line in enumerate(lines):
        if "[ dihedraltypes ]" in line and "improper" in line:
            i += 2
            while i < len(lines):
                impropertypes.append(tuple(lines[i].split()[:4]))
                i += 1
            break
    
    missing_dihedrals = []
    for dihedral in dihedrals:
        if dihedral[0] not in dihedraltypes and \
            dihedral[0][::-1] not in dihedraltypes:
            atom_field = " {:<4}{:<4}{:<4}{:<4}  ".format(*dihedral[0])
            params = "{}       {:>6}    {:>8}     {}".format(*dihedral[1])
            missing_dihedrals.append(atom_field + params)
            
    missing_impropers = []
    for improper in dihedrals_improper:
        if improper[0] not in impropertypes and \
            improper[0][::-1] not in impropertypes:
            atom_field = "{:<4}{:<4}{:<4}{:<4}     ".format(*improper[0])
            params = "{}      {:>6}    {:>8}     {}".format(*improper[1])
            missing_impropers.append(atom_field + params)  

    comment = "  ; added by Finn\n"    
    
    if missing_impropers != []:
        print("The improper dihedral(s): ")
        print(" i   j   k   l     func     phase      kd         pn")
        print(*missing_impropers, sep=comment, end=comment)
        print("are missing from the forcefield.")
        print("Add them to ffbonded under dihedral types.\n")
    
    if missing_dihedrals != []:
        print("The proper dihedral(s): ")
        print(" i   j   k   l     func     phase      kd         pn")
        print(*missing_dihedrals, sep=comment, end=comment)
        print("are missing from the forcefield.")
        print("Add them to ffbonded under dihedral types.\n")
        
def reorder_atoms():
    """
    Reorder atoms to follow style of Tyrosine in .rtp.
    
    Returns
    -------
    None.
    
    """       
    # dictionary of order for tyrosine-like atoms
    numbering = { "N" : 1, "H" : 2, "CA" : 3, "HA" : 4, "CB" : 5, "HB1" : 6, \
                 "HB2" : 7, "CG" : 8, "CD1" : 9, "HD1" : 10, "CE1" : 11, \
                 "HE1" : 12, "CZ" : 13, "OH" : 14, "HH" : 15, "CE2" : 16, \
                 "HE2" : 17, "CD2" : 18, "HD2" : 19, "C" : 20, "O" : 21 }
    mimetic = list(filter(lambda x : x.residue == residue_name, Atom.atoms))
    Atom.reordered = [None] * len(mimetic)
    non_tyr = []
    
    # reorder tyr-like atoms to the familiar index
    for mim in mimetic:
        if mim.name in numbering:
            Atom.reordered[numbering[mim.name] - 1] = mim
        else:
            non_tyr.append(mim)
            
    # add in the remaining atoms to reordered list vacancies
    i = 0
    for j, atom in enumerate(Atom.reordered):
        if atom == None:
            Atom.reordered[j] = non_tyr[i]
            i += 1
    
    for k, atom in enumerate(Atom.reordered):
        atom.num = k + 6
    
    return Atom.reordered

if __name__ ==  '__main__':
    main(argv)
    