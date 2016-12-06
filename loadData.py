# -*- coding: utf-8 -*-
"""
This script loads the SNP data from a pair of diseases and estimates the 
genetic correlation per chromosome

Proposed methodology:

For each chromosome:
    Region = Whole Chromosome
    

"""

import pandas as pd
import os
import numpy as np

rootdir = '../6.047-Data/'  # Change this to the path to the data when running on your machine
disease1_SNP = 'pgc.bip.full.2012-04.txt'
disease2_SNP = 'pgc.bip.full.2012-04.txt'

# Load data into a pandas df

disease1 = pd.read_csv(rootdir+disease1_SNP, sep='\t',names=['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'or', 'se', 'pval', 'info', 'ngt', 'CEUaf'],skiprows=[0]) 
disease2 = pd.read_csv(rootdir+disease2_SNP, sep='\t',names=['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'or', 'se', 'pval', 'info', 'ngt', 'CEUaf'],skiprows=[0]) 

# Convert the relevant columns to numeric values
disease1[['bp']] = disease1[['bp']].apply(pd.to_numeric, errors='coerce')
disease2[['bp']] = disease2[['bp']].apply(pd.to_numeric, errors='coerce')


print 'Loaded data'
print '\t Disease 1: \n'
print disease1.head()
print '\t Disease 2: \n'
print disease2.head()

chromosomes = set(disease1.hg18chr + disease2.hg18chr)

def estimate_corr(chromosome, region_start, region_end):
    '''Given a chromosome and region, creates a SNP data file for each disease,
    computes the genetic correlation between the two diseases in that region, 
    then deletes the created file.
    '''
    d1_file = disease1[disease1.hg18chr==chromosome][disease1.bp>=region_start][disease1.bp<=region_end]    
    d2_file = disease1[disease2.hg18chr==chromosome][disease2.bp>=region_start][disease2.bp<=region_end]

    # Save dataframes to textfile 
    filename = str(chromosome)+'_'+str(region_start)+'_'+str(region_end)
    d1_file.to_csv(rootdir+filename+'1.txt',sep='\t')
    d2_file.to_csv(rootdir+filename+'2.txt',sep='\t')
    
    # Estimate the genetic correlation
    corr = get_genetic_corr(rootdir+filename+'1.txt', rootdir+filename+'2.txt')
    
    # Remove files from folder
    os.remove(rootdir+filename+'1.txt')
    os.remove(rootdir+filename+'2.txt')
    
    return corr
    
def get_genetic_corr(disease1_file, disease2_file):
    '''Runs mungestat and ldsc on two diseases to estimate the genetic correlation'''
    # TODO: Complete this
    pass

def recursive_get_regions(chromosome, region_start, region_end):
    corr = estimate_corr(chromosome, region_start, region_end)
    # Base Case
    if corr < MIN_CORR or (region_end - region_start) <=1:
        return [(region_start, region_end, corr)]
    else:
        # Recurse and find the minimum correlation regions in the upper half and lower half
        region_mid = np.ceil((region_start+region_end)/2)
        return  recursive_get_regions(chromosome, region_start, region_mid) + recursive_get_regions(chromosome, region_mid, region_end)         

def get_minimal_regions(chromosome):
    '''Finds the minimal regions with significant genetic correlation'''
    
    # Initialize region start and rend
    region_start = 0
    region_end = 1
    
    # Find the minimal regions with sufficient genetic correlation
    minimal_regions = recursive_get_regions(chromosome, region_start, region_end)
    
    return minimal_regions

def compute_minimal_regions():
    ''' Computes the minimal regions with sufficient genetic correlation for 
    the two diseases'''
    
    minimal_regions = {}
    for chromosome in chromosomes:
        regions = get_minimal_regions(chromosome)
        minimal_regions[chromosome] = regions
    
    
    
