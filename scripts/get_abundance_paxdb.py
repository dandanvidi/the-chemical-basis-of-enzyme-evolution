import os
import pandas as pd

MapPath = '../data/paxdb-uniprot-links-v4_0/'
mapListing = sorted(os.listdir(MapPath))
paxdb2uniprot = pd.concat([pd.read_csv(MapPath+f, sep='\t', header=None, 
                           names=['paxdb_id', 'Entry']) for f in mapListing])

#%%
AbundancePath = '../data/paxdb-abundance-files-v4_0/'
abundanceListing = sorted(os.listdir(AbundancePath))
abundanceListing = {x:os.listdir(AbundancePath+x) for x in abundanceListing}
relevant_files = {}
for org_id, org_files in abundanceListing.iteritems():
    if len(org_files) == 1:
        relevant_files[org_id] = org_files[0]
    else:
        find = [i for i, f in enumerate(org_files) if 'WHOLE_ORGANISM-integrated.txt' in f]
        if find:
            relevant_files[org_id] = find[0]
            print org_files
