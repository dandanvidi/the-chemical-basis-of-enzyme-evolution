# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 16:34:10 2017

@author: dan
"""
import pandas as pd
import settings

iJO1366_abundance = settings.read_cache("iJO1366_2_abundance")
iJO1366_abundance = iJO1366_abundance[['reaction_id', 'glucose']]
iJO1366_abundance.columns = ['reaction_id', 'abundance[g/gCDW]']
iJO1366_kcat = settings.read_cache("iJO1366_2_kcat")
BRENDA_turnover = settings.read_cache("turnover")

median_kcat = BRENDA_turnover.groupby('EC_number', as_index=False).median()
#ORGANISM = 'Escherichia coli'

df = iJO1366_kcat.merge(median_kcat)
df = df.merge(iJO1366_abundance)

df['relative_rate'] = df['kcat [1/s]'] / df['Turnover_Number']

df.plot(kind='scatter', x='abundance[g/gCDW]', y='relative_rate', logx=True, logy=True)