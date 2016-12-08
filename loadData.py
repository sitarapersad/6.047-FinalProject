# -*- coding: utf-8 -*-
"""
This script loads the SNP data from a pair of diseases and estimates the 
genetic correlation per chromosome
"""

import pandas as pd
import os
import numpy as np
import subprocess
from StringIO import StringIO

rootdir = '../6.047-Data/'  # Change this to the path to the data when running on your machine
disease1_SNP = 'pgc.cross.BIP11.2013-05.txt'
disease2_SNP = 'pgc.cross.SCZ17.2013-05.txt'

global disease1
global disease2

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

DISEASES = {}
DISEASES['aut'] = {'sample_size': (4788 + 161, 4788 + 526), 'filename': 'pgc.cross.AUT8.2013-05.txt'}
DISEASES['add'] = {'sample_size': (1947 + 840, 1947 + 688), 'filename': 'pgc.cross.ADD4.2013-05.txt'}
DISEASES['bip'] = {'sample_size': (6990, 4820), 'filename': 'pgc.cross.BIP11.2013-05.txt'}
DISEASES['mdd'] = {'sample_size': (9227, 7383), 'filename': 'pgc.cross.MDD9.2013-05.txt'}
DISEASES['scz'] = {'sample_size': (9379, 7736), 'filename': 'pgc.cross.SCZ17.2013-05.txt'}

def munge(disease1, disease2):
    ''' Runs mungestat on two diseases in preparation for ldsc'''

    disease1_file = DISEASES[disease1]
    disease2_file = DISEASES[disease2]
    
    subprocess.call(['python', 'ldsc/munge_sumstats.py',
                     '--sumstats', '../6.047-Data/'+disease1_file,
                     '--N', '11810',
                     '--out', str(disease1),
                     '--merge-alleles', '../6.047-Data/w_hm3.snplist']
                    )

    subprocess.call(['python', 'ldsc/munge_sumstats.py',
                     '--sumstats', '../6.047-Data/'+disease2_file,
                     '--N', '17115',
                     '--out', str(disease2),
                     '--merge-alleles', '../6.047-Data/w_hm3.snplist']
                    )
    return None

def get_genetic_corr(disease1, disease2):
    '''Runs ldsc on two diseases to estimate the genetic correlation
    between the two selected diseases. Assumes that the data has already
    been munged and that the relevant .sumstats.gz files have been created'''

    # Run ldsc 
    subprocess.call(['python', 'ldsc/ldsc.py',
                     '--rg', str(disease1)+'.sumstats.gz',str(disease2)+'.sumstats.gz',
                     '--ref-ld-chr', '../6.047-Data/eur_w_ld_chr/',
                     '--w-ld-chr', '../6.047-Data/eur_w_ld_chr/',
                     '--out', str(disease1)+'_'+str(disease2)]
                    )
                    
    f = open(str(disease1)+'_'+str(disease2)+'.log', 'r')
    for line in f:
        if line =='Summary of Genetic Correlation Results\n':
            break
    total = 0
    lines = []
    for line in f:
        if total == 2:
            break
        lines.append(' '.join(s for s in line.split(' ') if s!= ''))
        total+=1

    f.close()
    summary = StringIO(''.join(lines))

    df = pd.read_csv(summary, sep=" ")

    to_remove = [str(disease1)+'_'+str(disease2)+'.log']
    for file in to_remove:
        os.remove(file)
    return float(df['rg'])

def get_SNPS_in_range(region_start, region_end, disease_df, chromosome):
    ''' Returns a set of SNP ids within the desired range on the given chromosome 
    in a disease dataframe.
    '''
    d1_file = disease_df[disease_df.hg18chr==chromosome][disease_df.bp>=region_start][disease_df.bp<=region_end]
    SNPs = d1_file['snpid']
    return set(SNPs)
    
    
def partition_sumstats(sumstats_df, SNP_set):
    ''' Given a sumstats_df and a set of target SNPids, return the 
    sumstats dataframe with the desired partition only.
    '''
    ###     # Load sumstats data into dataframe
    #sumstats_df = pd.read_table(sumstats_name+'.sumstats.gz', compression='gzip')
    
    # Get the rows that are in the region we want, based on the set of SNPs   
    return sumstats_df.loc[sumstats_df['SNP'].isin(SNP_set)]
    

def estimate_corr(chromosome, region_start, region_end):
    '''Given a chromosome and region, creates a SNP data file for each disease,
    computes the genetic correlation between the two diseases in that region, 
    then deletes the created file.
    '''
    
    
    # Estimate the genetic correlation
    corr = get_genetic_corr(filename+'1.txt', filename+'2.txt')
    
    # Remove files from folder
    os.remove(rootdir+filename+'1.txt')
    os.remove(rootdir+filename+'2.txt')
    return corr


def recursive_get_regions(chromosome, region_start, region_end):
    corr = estimate_corr(chromosome, region_start, region_end)
    # Base Case
    # TODO: How do we determine the MIN_CORR?
    if corr < MIN_CORR or (region_end - region_start) <=1:
        return [(region_start, region_end, corr)]
    else:
        # Recurse and find the minimum correlation regions in the upper half and lower half
        region_mid = np.ceil((region_start+region_end)/2)
        return  recursive_get_regions(chromosome, region_start, region_mid) + recursive_get_regions(chromosome, region_mid, region_end)         


def get_minimal_regions(chromosome):
    '''Finds the minimal regions with significant genetic correlation'''
    
    # Initialize region start and rend
    region_start = disease1[disease1.hg18chr==chromosome].bp.min()
    region_end = disease1[disease1.hg18chr==chromosome].bp.max(),
    
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
    


