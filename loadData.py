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

DISEASES = {}
DISEASES['aut'] = {'sample_size': (4788 + 161, 4788 + 526), 'filename': 'pgc.cross.AUT8.2013-05.txt'}
DISEASES['add'] = {'sample_size': (1947 + 840, 1947 + 688), 'filename': 'pgc.cross.ADD4.2013-05.txt'}
DISEASES['bip'] = {'sample_size': (6990, 4820), 'filename': 'pgc.cross.BIP11.2013-05.txt'}
DISEASES['mdd'] = {'sample_size': (9227, 7383), 'filename': 'pgc.cross.MDD9.2013-05.txt'}
DISEASES['scz'] = {'sample_size': (9379, 7736), 'filename': 'pgc.cross.SCZ17.2013-05.txt'}




global disease1_df
global disease2_df

# Load data into a pandas df
disease1_df = pd.read_csv(rootdir+disease1_SNP, sep='\t',names=['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'or', 'se', 'pval', 'info', 'ngt', 'CEUaf'],skiprows=[0]) 
disease2_df = pd.read_csv(rootdir+disease2_SNP, sep='\t',names=['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'or', 'se', 'pval', 'info', 'ngt', 'CEUaf'],skiprows=[0]) 

# Convert the relevant columns to numeric values
disease1_df[['bp']] = disease1_df[['bp']].apply(pd.to_numeric, errors='coerce')
disease2_df[['bp']] = disease2_df[['bp']].apply(pd.to_numeric, errors='coerce')

print 'Loaded data'
print '\t Disease 1: \n'
print disease1_df.head()
print '\t Disease 2: \n'
print disease2_df.head()

chromosomes = set(disease1_df.hg18chr + disease2_df.hg18chr)



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

def get_genetic_corr(disease1_sumstats_file, disease2_sumstats_file):
    '''Runs ldsc on two diseases to estimate the genetic correlation
    between the two selected diseases. Assumes that the data has already
    been munged and that the relevant .sumstats.gz files for the desired regions 
    have been created.'''

    # Run ldsc 
    subprocess.call(['python', 'ldsc/ldsc.py',
                     '--rg', disease1_sumstats_file+','+disease2_sumstats_file,
                     '--ref-ld-chr', '../6.047-Data/eur_w_ld_chr/',
                     '--w-ld-chr', '../6.047-Data/eur_w_ld_chr/',
                     '--out', disease1_sumstats_file+'_'disease2_sumstats_file]
                    )
                    
    f = open(disease1_sumstats_file+'_'disease2_sumstats_file+'.log', 'r')
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

    to_remove = [disease1_sumstats_file+'_'disease2_sumstats_file+'.log']
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
    '''Given a chromosome and region, creates a sumstats data file for each disease,
    computes the genetic correlation between the two diseases in that region, 
    then deletes the created file.
    '''
    
    # Get a list of the desired SNPids
    snpids1 = get_SNPS_in_range(region_start, region_end, disease1_df, chromosome)
    snpids2 = get_SNPS_in_range(region_start, region_end, disease1_df, chromosome)
    
    # Create the .sumstats.gz files for these chosen snpids
    partition1 = partition_sumstats(disease1_sumstats_df, snpids1)
    partition2 = partition_sumstats(disease2_sumstats_df, snpids2)
 
    # Save these to a .gz file
    sumstats_filename1 = '.sumstats.gz'
    sumstats_filename2 = '.sumstats.gz'
    partition1.to_csv(sumstats_filename1, '\t', compression='gzip')
    partition2.to_csv(sumstats_filename2, '\t', compression='gzip')
     
    # Estimate the genetic correlation
    corr = get_genetic_corr(sumstats_filename1, sumstats_filename2)
    
    # Remove the .gz files
    os.remove(rootdir+filename+'1.txt')
    os.remove(rootdir+filename+'2.txt')
    return corr


def recursive_get_regions(chromosome, region_start, region_end):
    ''' Given a chromosome and a start, stop coordinate, estimate the shared heritability between those two 
    diseases in that region.
    '''
    #Estimate the correlation
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
    disease1_region_start = disease1_df[disease1_df.hg18chr==chromosome].bp.min()
    disease1_region_end = disease1_df[disease1_df.hg18chr==chromosome].bp.max(),
    disease2_region_start = disease2_df[disease2_df.hg18chr==chromosome].bp.min()
    disease2_region_end = disease2_df[disease2_df.hg18chr==chromosome].bp.max(),
    
    region_start = min(disease1_region_start, disease2_region_start)
    region_end = max(disease1_region_end, disease2_region_end)
    
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
    


