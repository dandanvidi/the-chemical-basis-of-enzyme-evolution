# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 10:24:14 2017

@author: dan
"""
import pandas as pd
import settings
import numpy as np
from component_contribution.kegg_reaction import KeggReaction

manual_additions = settings.read_data('KEGG_manual_fixes')
KEGG_COMPOUNDS = settings.read_data('kegg2chebi')
KEGG_COMPOUNDS = pd.concat([KEGG_COMPOUNDS, manual_additions])
KEGG_COMPOUNDS = KEGG_COMPOUNDS['name'].to_dict()
KEGG_REACTIONS = pd.read_json("../data/kegg_reactions.json")
KEGG_REACTIONS['ECS'] = map(lambda x: "|".join(x), KEGG_REACTIONS['ECS'])
KEGG_REACTIONS = KEGG_REACTIONS.rename(columns={'ECS':'EC_number'})
KEGG_REACTIONS = settings.tidy_split(KEGG_REACTIONS, 'EC_number')

# manual EC fixes
KEGG_REACTIONS.EC_number[KEGG_REACTIONS.RID=='R07520'] = '1.14.99.27'
# Deleted entries
unspecified_ECS = ['-.-.-.-']
obsolete_enzymes = ['4.2.1.61']
remove_entries = unspecified_ECS+obsolete_enzymes
KEGG_REACTIONS = KEGG_REACTIONS[~KEGG_REACTIONS.EC_number.isin(remove_entries)]

KEGG_REACTIONS.set_index('EC_number', inplace=True)

sprs = ["|".join([";".join(map(str,x)) for x in tup]) for tup in KEGG_REACTIONS['reaction']]
reac_dict = [[list(reversed(x)) for x in tup] for tup in KEGG_REACTIONS['reaction']]
reac_dict = [dict(x) for x in reac_dict]
KEGG_REACTIONS['KEGG_string'] = reac_dict
KEGG_REACTIONS['human_readable'] = [{KEGG_COMPOUNDS[k]:v for k,v in d.iteritems()} for d in reac_dict]
KEGG_REACTIONS['human_readable'] = [KeggReaction(x).write_formula() for x in KEGG_REACTIONS.human_readable]
sprs = ["|".join([";".join(map(str,x)) for x in tup]) for tup in KEGG_REACTIONS['reaction']]
KEGG_REACTIONS['reaction'] = sprs
KEGG_REACTIONS = settings.tidy_split(KEGG_REACTIONS, 'reaction')
KEGG_REACTIONS[['stoichiometry', 'KEGG_ID']] = KEGG_REACTIONS['reaction'].str.split(';', expand=True)
KEGG_REACTIONS['stoichiometry'] = KEGG_REACTIONS['stoichiometry'].astype(float)

KEGG_REACTIONS['direction'] = KEGG_REACTIONS['stoichiometry'].apply(np.sign)

KEGG_REACTIONS['names'] = KEGG_REACTIONS['names'].apply(lambda x: ";".join(x))
KEGG_REACTIONS.drop('reaction', axis=1, inplace=True)

kegg_reverse = KEGG_REACTIONS.copy()
kegg_reverse.stoichiometry *= -1.
kegg_reverse.direction *= -1.
KEGG_REACTIONS = pd.concat([KEGG_REACTIONS, kegg_reverse])
KEGG_REACTIONS['EC_number'] = KEGG_REACTIONS.index
KEGG_REACTIONS.index = range(len(KEGG_REACTIONS))
KEGG_REACTIONS.to_csv('../cache/kegg_reactions.csv')
