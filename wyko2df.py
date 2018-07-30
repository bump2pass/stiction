# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 10:35:39 2018

@author: rob.macdonald
"""

import pandas as pd
from itertools import islice
import os
from os.path import  join

#list files in a directory
dr = 'C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\L18A183E'
lf = os.listdir(dr)
lf = [join(dr, i) for i in lf]
# a particular file to play with
fn = lf[0]
def wyko2df(fn):
    #read in the array containing the data
    df =  pd.read_csv(fn, delimiter = ',', header=None, skiprows =13, nrows = 736,na_values = ["Bad", "BAD", "bad"],
                        skip_blank_lines =  False, dtype = 'float64')
    ar = df.values#convert to array
    #read in the wavelength and multi paramters
    #first get all the header data, first 13 lines
    with open(fn) as myfile:
        head = list(islice(myfile, 13))
    #strip newline character from end of file
    head = [x.rstrip('\n') for x in head]
    #mult paramter is a fudge factor not explained in manual
    #other than you need to multiply by it get the correct height data
    #the following assumes that mult is always in the same place in the 
    #data structure. could be wrong!!!
    mult = head[12]
    mult = float(mult.split(',')[3])
    #similarly for wavelength
    wavelength = head[10]
    wavelength = float(wavelength.split(',')[3])
    #now multiply each value in the array ar
    #wavelength is nm/wave so we conver to um/wave by div 1000
    ar = ar*wavelength/mult/1000
    return ar




#read in each file in test directory
