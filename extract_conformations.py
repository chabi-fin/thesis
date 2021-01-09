#!/usr/bin/env python
# coding: utf-8

import numpy as np
import mdtraj as md
import os

home = os.getcwd()

fTyr = md.load_trr("fTyr_md.trr", top="fTyr_md.gro").remove_solvent()

# each frame is 100 ps
configs = [c for c, t in zip(fTyr, fTyr.time) if t % 1000 == 0]

for i, c in enumerate(configs):
    name = "config_" + str(i) 
    if not os.path.isdir(name):
        os.mkdir(name)
    file = name + "/geo.pdb"
    c.save_pdb(file)


os.chdir(home)



