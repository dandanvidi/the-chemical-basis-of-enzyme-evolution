# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 14:41:16 2017

@author: dan
"""
from bioservices.uniprot import UniProt
from io import StringIO
import time
import os
import pandas as pd

def get_enzyme_sequences(df):
    
    u = UniProt(verbose=False)
    lookups = set(df.groupby(['EC_number', 'Organism']).size().index)         
    report = pd.DataFrame(columns=['Entry','Sequence','Status','EC_number','Organism'])
    
    if os.path.isfile("../cache/enzyme_sequences.csv"):
        report = pd.DataFrame.from_csv("../data/reaction_to_genes_per_condition.pkl")
        done = set(report.groupby(['EC_number', 'Organism']).size().index)
        lookups = lookups.difference(done)
        

    for i, l in enumerate(lookups):
        EC, org = l
        try:
            params = ' '.join(['brenda:%s'%EC, 'organism:%s'%org]) 
            upids = u.search(params, columns="id, sequence, reviewed")
            x = pd.read_csv(StringIO(upids), sep='\t')
            x['EC_number'] = EC
            x['Organism']  = org
            report = pd.concant([report, x])
            report.to_csv("../cache/enzyme_suquences.csv")
        except:
            print "FAILED TO FIND EC %s in %s" %(EC, org)
        print "%i out of %i" %(i, len(lookups))
        time.sleep(0.01)             

    return report

if __name__ == '__main__':
    
    DF = pd.DataFrame.from_csv("../cache/BRENDA_by_KEGG.tsv", sep='\t')    
    get_enzyme_sequences(DF)