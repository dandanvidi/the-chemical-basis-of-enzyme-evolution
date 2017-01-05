#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 17:48:51 2017

@author: dan
"""

import pandas as pd
import seaborn as sns
import numpy as np

#df = pd.DataFrame.from_csv("../cache/BRENDA_by_KEGG.tsv", sep='\t')
df = pd.DataFrame.from_csv("../data/turnover.csv")

df = df[df.EC_number=="4.2.1.1"]
# CO2 hydration (also an hydrolysis reaction in a way) - Carbonic anhydrases
#CO2_hydration="CO2 + H2O ⇌ H2CO3"
#EC = "4.2.1.1"
#knon
#
#df[['EC1', 'EC2', 'EC3']] = df[['EC1', 'EC2', 'EC3']].astype(str) 
#df['EC_cat'] = df.EC1.str.cat(df.EC2, sep='.')
#EC_families = ['4.2', '3.1', '3.5', '3.4']
#df = df[df.EC_cat.isin(EC_families)]
#
#df.Turnover_Number = np.log10(df.Turnover_Number)
#sns.violinplot(x="EC_cat", y="Turnover_Number", data=df)