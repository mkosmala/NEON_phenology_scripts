#!/usr/bin/python

import csv
import sys
import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt

# -*- coding: utf-8 -*-
"""
Create figure for paper


Created on Wed May  9 13:57:27 2018

@author: mkosmala
"""



########
# MAIN #
########

if len(sys.argv) < 2:
    print ("format: create_variation_by_phenophase_figure.py <stats file>")
    exit(1)

infilename = sys.argv[1]

data = []
with open(infilename,'r') as infile:
    ireader = csv.reader(infile)
    
    # header
    next(ireader)
    
    # stats
    for row in ireader:
        data.append(row)
        
# plot
fig, ax = plt.subplots(1,1)

for trans in data:
    s = float(trans[3])
    scale = np.exp(float(trans[1]))
    x = np.linspace(lognorm.ppf(0.0001,s,0,scale),
                    lognorm.ppf(0.9999,s,0,scale), 10000)
    ax.plot(x, lognorm.pdf(x,s,0,scale), '-', alpha=0.6,
            label=trans[0])

frame1 = plt.gca()
frame1.axes.get_yaxis().set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.xlim(-2,42)
plt.xlabel("Standard deviation (days)")
plt.legend(frameon=False)

#plt.show()
plt.savefig('distribution_of_standard_deviations.png', bbox_inches='tight')


