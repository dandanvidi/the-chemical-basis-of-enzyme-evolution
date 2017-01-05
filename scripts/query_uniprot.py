# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 12:19:54 2016

@author: dan
"""
from bioservices.uniprot import UniProt
import pandas as pd
from io import StringIO
import time


BarEven = "../source/bi2002289_si_003.xls"
kinetictable = pd.read_excel(BarEven, 1)
kcat = pd.read_csv("../data/turnover.csv").dropna(subset=['KEGG_ID'])
kM   = pd.read_csv("../data/km.csv").dropna(subset=['KEGG_ID'])
kcat_lookups = set(kcat.groupby(['EC_number', 'Organism']).size().index)
KM_lookups   = set(kM.groupby(['EC_number', 'Organism']).size().index)
lookups = kcat_lookups.union(KM_lookups)
u = UniProt(verbose=False)

##already_have_it = pd.DataFrame.from_csv("../cache/sequences.csv")
##already_couldnt_find = pd.DataFrame.from_csv("../cache/failed.csv")
##x = set(already_have_it.groupby(['EC_number', 'Organism']).size().index)
##y = set(already_couldnt_find.groupby(['EC_number', 'Organism']).size().index)
#
##not_look = x.union(y)
##lookups = lookups.difference(not_look)
#
#df = pd.DataFrame(columns=['Entry','Sequence','Status','EC_number','Organism'])
#failed = pd.DataFrame(columns=['EC_number','Organism'])
#for i, l in enumerate(lookups):
#    EC, org = l
#    try:
#        params = ' '.join(['brenda:%s'%EC, 'organism:%s'%org]) 
#        upids = u.search(params, columns="id, sequence, reviewed")
#        x = pd.read_csv(StringIO(upids), sep='\t')
#        x['EC_number'] = EC
#        x['Organism']  = org        
#        df = pd.concat([df, x])
#        df.to_csv("../cache/sequences.csv")
#    except ValueError:
#        tmp = pd.DataFrame.from_records([l], columns = ['EC_number', 'Organism'])
#        failed = pd.concat([failed,tmp])
#        failed = pd.DataFrame.from_records(failed)
#        failed.to_csv("../cache/failed.csv")
#    print "%i out of %i" %(i, len(lookups)) 
#    time.sleep(0.01)        
#
