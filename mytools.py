# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 08:09:57 2021

@author: Lauren
"""

from sys import argv, exit

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