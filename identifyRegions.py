# -*- coding: utf-8 -*-
"""
Main code for identifying regions of interest
"""

import loadData as ld
import pandas as pd
import numpy as np
import geneAnnotate as ga
import matplotlib.pyplot as plt

def find_regions(disease1, disease2, cutoff_pval=0.05):
    df = pd.read_csv('RESULTS_' + disease1 + '_' + disease2 + '.txt',
                     sep='\t', names=['chromosome', 'region_start', 'region_end', 'correlation', 'pvalue', 'std_err', 'SNPs']
                     )
    df = df.apply(pd.to_numeric, errors='coerce')

    layer = [0]
    for i in xrange(len(df)-1):
        row_1 = df.loc[i]
        row_2 = df.loc[i+1]
        if row_1['chromosome'] != row_2['chromosome']:
            layer.append(0)
        else:
            distance_1 = row_1['region_end'] - row_1['region_start']
            distance_2 = row_2['region_end'] - row_2['region_start']
            ratio = distance_2*1.0/distance_1
            layer.append(layer[-1] - int(np.round(np.log2(ratio))))
    df['layer'] = pd.Series(layer)
    correction = {}
    for layer in set(df['layer']):
        correction[layer] = len(df[df.layer==layer])

    reasonable = df[np.abs(df.correlation) < 1]
    reasonable = reasonable[reasonable.pvalue > 0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    offset = {1: 0}
    for i in xrange(2, 23):
        offset[i] = offset[i-1] + max(df[df.chromosome == (i-1)]['region_end'])
    for key, row in reasonable.iterrows():
        color = 0
        if (int(row['chromosome']) % 2) == 0:
            color = 0.5
        ax.plot([offset[int(row['chromosome'])]+row['region_start'],
                 offset[int(row['chromosome'])]+row['region_end']],
                [-np.log10(row['pvalue']), -np.log10(row['pvalue'])],
                '-', color=str(color), linewidth=1.5)
    ax.axhline(y=-np.log10(cutoff_pval/len(df)), xmin=0, xmax=1, color='r', linewidth=3, linestyle='--', markeredgecolor='k')
    plt.xticks([min(reasonable[reasonable.chromosome==1]['region_start'])+(offset[i-1]+offset[i])/2 for i in xrange(2, 23)]
               +[max(reasonable[reasonable.chromosome==22]['region_end'])/2+offset[22]], np.arange(1, 23),
                rotation='vertical')

    print reasonable.sort_values(['pvalue'])
    print cutoff_pval/len(df)
    plt.xlabel('Chromosome', fontsize=25)
    plt.ylabel('$-\log P$', fontsize=25)
    plt.tick_params(axis='both', labelsize=10)
    plt.ylim([0, -np.log10(min(reasonable['pvalue']))*1.1])
    plt.savefig('manhattan_'+disease1+'_'+disease2+'.pdf', bbox_inches='tight')
    output = reasonable[reasonable.pvalue <= cutoff_pval/len(df)]
    #return output.values.tolist()


if __name__ == '__main__':
    diseases = ['aut', 'add', 'bip', 'mdd', 'scz']
    for i in xrange(len(diseases) - 1):
        for j in xrange(i + 1, len(diseases)):
            print diseases[i], diseases[j]
            roi = find_regions(diseases[i], diseases[j])
            ga.region_to_bed('roi_bed/roi_'+diseases[i]+'_'+diseases[j]+'.bed', roi)