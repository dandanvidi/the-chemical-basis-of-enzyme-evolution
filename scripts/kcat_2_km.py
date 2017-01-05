import pandas as pd
from math import ceil
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from bioservices.uniprot import UniProt
from io import StringIO
import time


def despine(ax, fontsize=15):
    ax.tick_params(right=0, top=0, direction='out', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(ax.get_xlabel(), size=fontsize)
    ax.set_ylabel(ax.get_ylabel(), size=fontsize)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    
if __name__ == '__main__':
    
    u = UniProt(verbose=False)
    # load data from Bar-Even et al 2011
    kinetictable   = pd.DataFrame("../cache/BRENDA_by_KEGG.csv", sep='\t')
    ecs = list(kinetictable.columns[0:4])
    kinetictable['EC'] = kinetictable[ecs].astype(str).apply(lambda x: '.'.join(x), axis=1)
    kcat = kinetictable.drop(['KM', 'Compound Id (KEGG)'], 1)
    km   = kinetictable.drop('kcat', 1)
    
    kcat.dropna(subset=['kcat'], inplace=True)
    kcat = kcat.groupby(ecs+['EC', 'Organism', 'Reaction direction (KEGG)',
                         'pH', 'Temperature'], as_index=False).median()

    km   = km.groupby(['EC', 'Organism', 'Compound Id (KEGG)',
                         'pH', 'Temperature'], as_index=False).median()

    df = km.merge(kcat, on=ecs+['EC', 'Organism', 'Reaction direction (KEGG)',
                            'pH', 'Temperature'], how='outer', sort=True)

    # group by EC and sort by amount of data
    # EC decimal sets the granulation depth    
    df.dropna(inplace=True)

'''
    X = 5
    EC_decimal = 4
    groups = df.groupby(ecs[0:EC_decimal], as_index=False)
    freqs = groups.size().sort_values(ascending=False)
    sorted_groups = [groups.get_group(i) for i,j in freqs.iteritems()]    
    sdf = pd.concat(i for i in sorted_groups[:X])
    
    lookups = set(sdf.groupby(['EC', 'Organism']).size().index)         


    failed = pd.DataFrame(columns=['EC','Organism'])
    report = []
    for i, l in enumerate(lookups):
        EC, org = l
        try:
            params = ' '.join(['brenda:%s'%EC, 'organism:%s'%org]) 
            upids = u.search(params, columns="id, sequence, reviewed")
            x = pd.read_csv(StringIO(upids), sep='\t')
            x['EC'] = EC
            x['Organism']  = org
            report.append(x)
#            df = df.merge(x, 'left')
#            df.to_csv("../cache/kinetic_table.csv")
        except ValueError:
            tmp = pd.DataFrame.from_records([l], columns = ['EC_number', 'Organism'])
            failed = pd.concat([failed,tmp])
            failed = pd.DataFrame.from_records(failed)
            failed.to_csv("../cache/failed.csv")
            print "FAILED TO FIND EC %s in %s" %(EC, org)
        
        print "%i out of %i" %(i, len(lookups))              
        time.sleep(0.01)        
    
    report = pd.concat(report)
    report = report[report['Status']=='reviewed']
    df = df.merge(report)
    df.to_csv("../cache/kinetic_table.csv")
#    # merge measurements from different publications - take the median
#    # use EC number, substrate and organism to identify the reaction
#    df = kinetictable.dropna(subset=kinetictable.columns[-2:], how='all')
#    merged_col = [u'EC1', u'EC2', u'EC3', u'EC4', u'Compound Id (KEGG)',
#       u'Reaction ID (KEGG)', u'Reaction direction (KEGG)', 'Organism ID','Temperature', 'pH']
#    median = df.groupby(merged_col, as_index=False).median()
#    median = median[merged_col+['kcat', 'KM']]
#    median.loc[:, 'log10_kcat'] = np.log10(median['kcat'])
#    median.loc[:, 'log10_KM'] = np.log10(median['KM'])
#    

#
#    sdf['EC'] = sdf[ecs].astype(str).apply(lambda x: '.'.join(x), axis=1)
#    sdf.replace(to_replace={'Organism ID':organismmaper}, inplace=True)
#    fill = '.-'*(4-len(ecs))
#    sdf['EC'] = sdf['EC'].apply(lambda x: x+fill)
#
#    lookups = set(sdf.groupby(['EC', 'Organism ID']).size().index)

#

#    # use only reactions that have kcat or km mneasured
#    
#    # filter extreme conditions - pH and temperature
#    df = df[(pd.isnull(df['Temperature'])) | 
#            (df['Temperature']<=37) & 
#            (df['Temperature']>=25)]
#    df = df[(pd.isnull(df['pH'])) | 
#            (df['pH']<=8) & 
#            (df['pH']>=6)]
#
#
#    g = sns.lmplot('log10_kcat', 'log10_KM', sdf, row='EC', fit_reg=False,
#                   hue='Compound Id (KEGG)', legend=False)
#    g = (g.set_axis_labels(r"$\log _{10} \ k_{cat} \left[s^{-1} \right]$", 
#                           r"$\log _{10} \ K_{M} \left[\mu M \right]$")
#        .set(xlim=(-3, 6), ylim=(-2, 6)))
#
#    g.savefig("../res/kcat2km_EC_decimal=%i.pdf"%len(ecs))

'''