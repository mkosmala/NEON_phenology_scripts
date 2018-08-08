#!/usr/bin/python

import csv
import sys
import matplotlib.pyplot as plt
import pandas as pd

# -*- coding: utf-8 -*-
"""
Visualization for output of how_many_individuals.py

Created on Mon May  7 13:26:58 2018

@author: mkosmala
"""

quants = [0.25,0.5,0.75,0.95,1.0]


def listtostr(inlist):
    return '-'.join(inlist)

# straight conversion
def xylist_from_dict(adict):
    xlist = []
    ylist = []
    for x in sorted(adict):
        xlist.append(x)
        ylist.append(adict[x])
    return xlist, ylist

# conversion with normalization
def normalize_from_dict(adict):
    normalizer = adict[30]
    xlist, ylist = xylist_from_dict(adict)
    ynorm = [y / normalizer for y in ylist]
    return xlist, ynorm

# normalize = T/F, data = mean or stdev
def tabulate(normalize,data):
    p5 = []
    p10 = []
    p15 = []
    p20 = []
    p25 = []
    p30 = []
    for phenoset in data:
        d = data[phenoset]
        if normalize:            
            p5.append(d[5]/d[30])
            p10.append(d[10]/d[30])
            p15.append(d[15]/d[30])
            p20.append(d[20]/d[30])
            p25.append(d[25]/d[30])
            p30.append(d[30]/d[30])
        else:
            p5.append(d[5])
            p10.append(d[10])
            p15.append(d[15])
            p20.append(d[20])
            p25.append(d[25])
            p30.append(d[30])
        
    # take the quantiles
    df = pd.DataFrame()
    for k,name in zip([p5,p10,p15,p20,p25,p30],
                      ['5','10','15','20','25','30']):
        ser = pd.Series(k).quantile(quants)
        colnames = df.columns
        df = df.assign(q=ser.values)
        df.columns = list(colnames) + [name]

    # name the rows
    df = df.assign(quantiles=quants)
    df = df.set_index('quantiles')
 
    return df


########
# MAIN #
########

if len(sys.argv) < 2:
    print ("format: viz_how_many_individuals.py <data file>")
    exit(1)

inputfilename = sys.argv[1]

# read in the data
mean_data = {}
std_data = {}
with open(inputfilename,'r') as infile:
    ireader = csv.reader(infile)
    
    # header
    next(ireader)
    
    for row in ireader:
        pheno, site, species, year, k, smean, sstd = row
        pkey = listtostr([pheno,site,species,year])
        if pkey not in mean_data:
            mean_data[pkey] = {}
        if pkey not in std_data:
            std_data[pkey] = {}
        mean_data[pkey][int(k)] = float(smean)
        std_data[pkey][int(k)] = float(sstd)
        

# quick stats for a table
# we want to know:
#    for 5, 10, 15, 20, 25, 30 individuals
#    - what is the max (normalized and unnormalized) lost precision
#    - for 25%, 50%, 75%, 95%, 100% of the data
#    - for both mean and std
# that's 4 tables
tab = tabulate(True,mean_data)
tab.to_csv("table_of_normalized_means.csv")
tab = tabulate(False,mean_data)
tab.to_csv("table_of_unnormalized_means.csv")
tab = tabulate(True,std_data)
tab.to_csv("table_of_normalized_std_devs.csv")
tab = tabulate(False,std_data)
tab.to_csv("table_of_unnormalized_std_devs.csv")


# mean plot
plt.style.use('grayscale')
fig, ax = plt.subplots()
for phenoset in mean_data:
    x, y = normalize_from_dict(mean_data[phenoset])
    plt.plot(x,y)

plt.xlabel("Number of individuals")
plt.ylabel("Increase in standard devation")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#plt.show()
plt.savefig('change_in_precision_of_mean.png', bbox_inches='tight')


# std plot
plt.clf()
for phenoset in std_data:
    x, y = normalize_from_dict(std_data[phenoset])
    plt.plot(x,y)

plt.xlabel("Number of individuals")
plt.ylabel("Increase in standard devation")

#plt.show()
plt.savefig('change_in_precision_of_std.png', bbox_inches='tight')



