# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 00:01:12 2016

@author: sitara
"""

"""
Code for selecting genes from regions.
"""

def region_to_bed(bed_name, regions):
    ''' Given a list of regions of the form [(chromosome,start,stop)], returns all the chromosomal positions in the region'''
    bedout = open(bed_name, 'w+')
    for region in regions:
        chromosome,start,stop = region
        chromosome = int(chromosome)
        start = int(start)
        stop = int(stop)
        bedout.write('chr'+str(chromosome)+'\t'+str(start)+'\t'+str(stop)+'\n')
    bedout.close()
    return