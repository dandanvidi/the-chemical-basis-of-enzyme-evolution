# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 14:45:18 2017

@author: dan
"""

import pandas as pd
import settings

def merge_KEGG_and_BRENDA(BRENDA_DF):

    kegg = settings.read_cache('kegg_reactions')    
    # drop kientic parameters with negative values    
    BRENDA_DF = BRENDA_DF[BRENDA_DF>0]    

    # drop mutated enzymes
    BRENDA_DF = BRENDA_DF[(pd.isnull(BRENDA_DF['Commentary'])) |
            (BRENDA_DF['Commentary'].str.find('muta') == -1)]
            
    kegg = kegg[kegg.stoichiometry<0]

    # keep biologically relevant reactions by inner merge BRENDA with KEGG
    DF = BRENDA_DF.merge(kegg, how='inner')

    # split EC to allow filterring per EC identifier
    
    return DF
    
        
if __name__ == '__main__': 



    turnover = pd.DataFrame.from_csv("../data/turnover.csv")
    turnover = merge_KEGG_and_BRENDA(turnover)    
    turnover = turnover.groupby(['EC_number', 'RID', 
                                 'Organism', 'direction'], as_index=False).median()

    turnover.drop(['stoichiometry', 'direction', 'LigandID'], axis=1, inplace=True)
    turnover.to_csv("../cache/turnover.csv")
#    km = pd.DataFrame.from_csv("../data/km.csv")
    
#    km = merge_KEGG_and_BRENDA(km, kegg)
#    km = km.groupby(['EC1', 'EC2', 'EC3', 'EC4', 'EC_number', 'RID', 
#                     'Organism', 'LigandID'], as_index=False).median()

#%%

#    turnover.drop(['stoichiometry', 'LigandID'], axis=1, inplace=True)

#    DF = km.merge(turnover, how='outer')
#    DF.to_csv("../cache/BRENDA_by_KEGG.tsv", sep='\t')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    