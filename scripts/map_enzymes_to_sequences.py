import pandas as pd
import numpy as np


def merge_KEGG_and_BRENDA(BRENDA_DF, KEGG_REACTIONS):
    
    # drop kientic parameters with negative values    
    BRENDA_DF = BRENDA_DF[BRENDA_DF>0]    

    # drop mutated enzymes
    BRENDA_DF = BRENDA_DF[(pd.isnull(BRENDA_DF['Commentary'])) |
            (BRENDA_DF['Commentary'].str.find('muta') == -1)]
            
    # keep biologically relevant reactions by inner merge BRENDA with KEGG
    DF = BRENDA_DF.merge(KEGG_REACTIONS)

    # split EC to allow filterring per EC identifier
    DF = split_ECnumber(DF)
    
    return DF
    
def split_ECnumber(df):
    df[['EC1', 'EC2', 'EC3', 'EC4']] = df['EC_number'].str.split(".", expand=True)
    return df
        
def sort_by_EC_frequency(df, EC_decimal=4):
    groups = df.groupby(['EC%s' %(i+1) for i in range(EC_decimal)])
    freqs = groups.size().sort_values(ascending=False)
    size_sorted_groups = [(i,j,groups.get_group(i)) for i, j in freqs.iteritems()]
    df = pd.concat(i[2] for i in size_sorted_groups)
    return df, freqs, size_sorted_groups    

if __name__ == '__main__': 


    kegg = pd.DataFrame.from_csv("../cache/kegg_reactions.csv", sep='\t')
    turnover = pd.DataFrame.from_csv("../data/turnover.csv")
    km = pd.DataFrame.from_csv("../data/km.csv")
    
    km = merge_KEGG_and_BRENDA(km, kegg)
    km = km.groupby(['EC1', 'EC2', 'EC3', 'EC4', 'EC_number', 'RID', 
                     'Organism', 'LigandID'], as_index=False).median()

    turnover = merge_KEGG_and_BRENDA(turnover, kegg)    
    turnover = turnover.groupby(['EC1', 'EC2', 'EC3', 'EC4', 'EC_number', 'RID', 
                                 'Organism', 'direction'], as_index=False).median()
    turnover.drop(['stoichiometry', 'LigandID'], axis=1, inplace=True)

    DF = km.merge(turnover, how='outer')
    DF.to_csv("../cache/BRENDA_by_KEGG.tsv", sep='\t')


#    df, freqs, size_sorted_groups = sort_by_EC_frequency(df, EC_decimal=4)
#    map(str,x)
#    l = []
#    for i in np.arange(1,6):
#        seqs = pd.read_csv("~/Documents/cache/%i__BRENDA_sequences.csv"%i, sep='\t', 
#                             skiprows=[0], dtype=str)
#        seqs.set_index('Accession_Code', inplace=True)
#    
#        l.append(df.merge(seqs))
#
#DF = pd.concat(l)
#DF.to_csv("../cache/complete_enzyme_table.csv")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    