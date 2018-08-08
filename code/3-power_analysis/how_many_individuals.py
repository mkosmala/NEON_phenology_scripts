#!/usr/bin/python

import csv
import sys
import numpy as np
from random import choices

# -*- coding: utf-8 -*-
"""
"Part 5"
How many individuals do we really need?
Uses only site-spp-year-phenophases that have 25+ individuals

Created on Mon May  7 11:31:12 2018

@author: mkosmala
"""

# Take a list of doys and for all subsets of size 5 to 30,
# bootstrap-draw a list of 10,000 doys with replacement.
# For each draw, calculate the mean and standard deviation
# returns two dicts keyed by subset size
def analyze(datelist):
    
    mean_data = {}
    std_data = {}
    
    # draw sets of 5 to 30 doys
    for k in range(5,31):
        
        # bootstrap replicates
        reps = []
        for i in range(0,1000):
            reps.append(choices(datelist,k=k))

        # calculate mean and stdev        
        draws = np.array(reps)
        means = np.mean(draws, axis=1)
        stds = np.std(draws, axis=1)

        # get the standard deviation of the means
        # and the standard deviation of the stdevs
        std_of_mean = np.std(means)
        std_of_std = np.std(stds)        

        mean_data[k] = std_of_mean
        std_data[k] = std_of_std
        
    return mean_data, std_data        

# check to make sure they're not all the same (i.e. SOME standard deviation)
def check_dates(datelist):
    if np.std(datelist) == 0:
        return False
    return True
    

def listtostr(inlist):
    return '-'.join(inlist)


########
# MAIN #
########

if len(sys.argv) < 4:
    print ("format: how_many_individuals.py <which combos file> <phenology file> <output file>")
    exit(1)

combosfilename = sys.argv[1]
phenofilename = sys.argv[2]
outputfilename = sys.argv[3]

# read in the combos file
combos = {}
with open(combosfilename,'r') as cfile:
    creader = csv.reader(cfile)

    # header    
    next(creader)
    
    for row in creader:
        site, spp, yr, gf, pp, n = row
        combos[listtostr([pp,site,spp,yr])] = [pp,site,spp,yr]
    
# read in the phenology file
# save in a data structure
# phenophase, site, spp, yr: list of dates as DOY
pdata = {}
with open(phenofilename,'r') as pfile:
    preader = csv.reader(pfile)
    
    # header
    next(preader)
    
    for row in preader:
        site = row[1]
        species = row[2]
        pheno = row[5]
        year = row[6]
        doy = int(row[7])
        
        pkey = listtostr([pheno,site,species,year])

        # only do phenology dates with enough individuals        
        if pkey in combos:

            if pkey not in pdata:
                pdata[pkey] = []
            pdata[pkey].append(doy)

# now the analyses
# track the standard deviations of the means and standard deviations of the
# bootstrap draws for each pkey
all_means = {}
all_stds = {}
for phenoset in pdata:
    dates = pdata[phenoset]
    if check_dates(dates):
        all_means[phenoset], all_stds[phenoset] = analyze(dates)

# output the data
with open(outputfilename,'w') as ofile:
    owriter = csv.writer(ofile)
    
    # header
    owriter.writerow(["phenophase","site","species","year","subset_size","std_of_mean","std_of_std"])
    
    # data
    for phenoset in all_means:
        means = all_means[phenoset]
        stds = all_stds[phenoset]
        for item in sorted(means):
            owriter.writerow(combos[phenoset] + [item,means[item],stds[item]])

