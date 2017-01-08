# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 13:37:50 2016

@author: noore

calculates dG'0 and dG'm for reaction from a Cobra model

"""

import settings
from ast import literal_eval
import numpy as np
import pandas as pd
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution_trainer \
    import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T

class reaction_thermodynamics(object):

    def __init__(self):

        reactions = set(settings.read_cache('kegg_reactions').KEGG_string.values)
        
        reactions = map(literal_eval, reactions)
        self.reactions = []
        self._not_balanced = []

        self.rstrings = []
        for r in reactions:
            k = KeggReaction(r)
            if k.is_balanced() and not k.is_empty():
                self.rstrings.append(k.write_formula())
                self.reactions.append(r)
            else:
                self._not_balanced.append(k)
        self.Kmodel = KeggModel.from_formulas(self.rstrings)

        self.cc = ComponentContribution.init()
        self.pH = 7.3
        self.I = 0.25
        self.RT = R * default_T

    def get_thermodynamics(self):
        '''
            Calculates the dG0 of a list of a reaction.
            Uses the component-contribution package (Noor et al) to estimate
            the standard Gibbs Free Energy of reactions based on
            component contribution approach and measured values
            (NIST and Alberty)

            Calculates the reversibility index (RI) of a reaction.
            The RI represent the change in concentrations of metabolites
            (from equal reaction reactants) that will make the reaction
            reversible. That is, the higher RI is, the more irreversible
            the reaction. A convenient threshold for reversibility is
            RI>=1000, i.e. a change of 1000-fold in metabolite concentrations
            is required in ordeer to flip the reaction direction.
        '''

        self.Kmodel.add_thermo(self.cc)

        temp = self.Kmodel.get_transformed_dG0(pH=self.pH, I=self.I, T=298.15)

        dG0_prime, dG0_cov = np.matrix(temp[0]), np.matrix(temp[1])

        conc = np.matrix(np.ones((len(self.Kmodel.cids), 1))) * 1e-3  # in M
        adj_mets = np.matrix(np.ones((len(self.Kmodel.cids), 1)))
        if 'C00001' in self.Kmodel.cids:
            j = self.Kmodel.cids.index('C00001')
            conc[j, 0] = 1
            adj_mets[j, 0] = 0

        dGm_prime = dG0_prime + self.RT * (self.Kmodel.S.T * np.log(conc))
        adj_count = np.abs(self.Kmodel.S.T) * adj_mets
        inv_adj_count = np.matrix(np.diag(
            map(lambda x: 1.0/x, adj_count.flat)))

        logRI = 2.0 * (inv_adj_count * dGm_prime) / (self.RT)

        res_df = pd.DataFrame(index=self.reactions)
        res_df['dG0_prime'] = dG0_prime
        res_df['dG0_prime_std'] = dG0_cov
        res_df['dGm_prime'] = dGm_prime
        res_df['logRI'] = logRI
        res_df['formula'] = self.rstrings
        res_df = res_df.round(2)
        return res_df

if __name__ == '__main__':
    Th = reaction_thermodynamics()
    df = Th.get_thermodynamics()
