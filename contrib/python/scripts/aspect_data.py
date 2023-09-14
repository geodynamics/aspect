#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Routines to read output files of ASPECT """

import numpy as np
import pandas as pd
import os

__author__ = 'The authors of the ASPECT code'
__copyright__ = 'Copyright 2023, ASPECT'
__license__ = 'GNU GPL 2 or later'



def read_statistics(fname):
    """ Read the statistics file output by CIG-ASPECT
    
    return a pandas table, where names are taken from the statistics file.
    """
    # header:
    header = []
    header_read = True

    with open(fname) as f:
        while header_read :
            line = f.readline()
            if line[0] == '#':
                idx_start = line.find(":")
                header.append(line[idx_start+2:-1])
            else:
                header_read = False
                
    # data
    values = pd.read_csv(fname, skiprows=len(header), header=None, delim_whitespace=True, names=header)
    return values



def read_gnuplot_visu(fname):
    """ Read a gnuplot file output by ASPECT for visualization 
    
    ! Only tested for single thread, 2D
    """
    # header:
    header = []
    header_read = True

    # Look for the line that starts with "# <", which contains all the names of variables. 
    with open(fname) as f:
        i = 0
        while header_read and i<100:
            line = f.readline()
            i += 1
            if line[0:3] == '# <':
                line_shorten = line[1:]
                idx_start = line_shorten.find("<")
                idx_end = line_shorten.find(">")
                i=0
                while idx_start >= 0 and  i< 100:
                    i+=1
                    header.append(line_shorten[idx_start+1:idx_end])
                    line_shorten = line_shorten[idx_end+1:]
                    idx_start = line_shorten.find("<")
                    idx_end = line_shorten.find(">")       
                header_read = False
    # search the duplicated names to add the suffixes _x _y.They can only be following each others. 
    old_name = [""]
    for i, name in enumerate(header): 
        if name == old_name: 
            header[i-1] = name+"_x"
            header[i] = name+"_y"
        old_name = name
        
    data = pd.read_csv(fname,comment='#', header=None, names=header, sep=" ", index_col=False)
    data.dropna(axis=0, how='any', thresh=None, subset=None, inplace=False) # remove the empty lines
    return data
