# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 09:20:20 2016

@author: dan
"""

from Bio.pairwise2 import align, _align, format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Only reactions with both kcat and KM
BRENDA = pd.DataFrame.from_csv("../cache/BRENDA_by_KEGG.tsv", sep='\t')
BRENDA.dropna(subset=['Turnover_Number', 'KM_Value'], inplace=True)
BRENDA['kcat/km'] = BRENDA['Turnover_Number'] / BRENDA['KM_Value']

#%%
# Merge enzyme sequences and kinetic data
df = pd.DataFrame.from_csv("../cache/enzyme_sequences.csv")
df = df[df.Status=='reviewed']
df = df.merge(BRENDA)

# analyze a specific EC class and group by reaction and metabolite
# for comparison
df = df[(df.EC1==2) & (df.EC2==7) & (df.EC3==1)]
DF = df.groupby(['RID', 'LigandID'])
#%%

# Calculate distance matrix for all enzymes that can catalyze the same reaction
# and have km and kcat data for a particular substrate

l = []
for it, (k, g) in enumerate(DF):
    g = g.drop_duplicates(subset=['Entry'])
    matrix = np.zeros((len(g), len(g)))
    for i, s1 in enumerate(g.Sequence):
        for j, s2 in enumerate(g.Sequence):
            if i<j:
                matrix[i,j] = align.globalms(s1, s2, 0, -1, -1, -1, score_only=True)
            elif i>j: 
                matrix[i, j] = matrix[j, i]
    index = g.Entry.copy()
    index.name = 'Entry_2'
    distance_df = pd.DataFrame(columns=g.Entry, index=index, data=matrix)
    distance_df = distance_df.unstack()
    distance_df.name = 'distance'
    distance_df = distance_df[distance_df<=0].reset_index()
    
    tmp = g.set_index('Entry')
    distance_df['Organism_2'] = tmp.loc[distance_df.Entry_2]['Organism'].values
    tmp = distance_df.groupby(['Organism_2', 'Entry'], as_index=False).max()
    k1 = g.set_index('Entry').loc[tmp.Entry, 'kcat/km']
    k2 = g.set_index('Entry').loc[tmp.Entry_2, 'kcat/km']
    y = k1.values / k2.values
    x = tmp.distance.values
    xy = pd.DataFrame(columns=['EC_number', 'Organism_1', 'Organism_2', 
                               '$\Delta k$', 'sequence distance'],
                      index=range(len(x)))
    xy['sequence distance'] = x
    xy['$\Delta k$'] = y
    xy['EC_number'] = g.EC_number.values[0]
    xy['Organism_1'] = g.Organism.values[0]
    xy['Organism_2'] = tmp.Organism_2.values[0]
    xy['Entry_1'] = k1.index
    xy['Entry_2'] = k2.index
    
    l.append(xy)

l = pd.concat(l)   
#%%

#f = sns.lmplot('$\Delta k$', 'sequence distance', l, fit_reg=False,
#               hue='EC_number', legend=False, logx=True)

l.plot(kind='scatter', x='sequence distance', y='$\Delta k$', logy=True)
plt.savefig("../delta k and distance.pdf")
#    if it > 5: break