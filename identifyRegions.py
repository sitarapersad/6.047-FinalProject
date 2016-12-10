# -*- coding: utf-8 -*-
"""
Main code for identifying regions of interest
"""

import loadData as ld
import pandas as pd
import numpy as np

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

    print df.head()
    print correction
    return df[df.pvalue <= pd.Series([cutoff_pval/correction[layer] for layer in df.layer])][df.pvalue > 0][np.abs(df.correlation) < 1].sort_values(by=['chromosome', 'region_start'])

