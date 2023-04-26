import sys

#sys.path.append('dspyields/kegg_helper')

import os
import json
from pathlib import Path
import cobra
# import cobra.test
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite
import pubchempy as pcp
from dspyields.kegg_helper.pykegg import *
from dspyields.utils import get_config
from dspyields.utils import read_json

cfg = get_config()

# if os.path.exists(dir + 'DatabaseVersions.csv'):
#     os.remove(dir + 'DatabaseVersions.csv')
# if os.path.exists('/path/to/directory/DatabaseVersions.csv')

rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])


def add_cpd2rxn(model, rxn):
    r = rxn_dict[rxn]
    rxn_coff = parser_equation_coef(r['equation'])
    reactants, products = split_equation(r['equation'])

    names = locals()

    meta_dict = {}
    for cpd in reactants:
        cpd = re.search(r'.*(C\d{5}).*', cpd).group(1)
        # c = MyCompound(cpd)
        c = cpd_dict[cpd]

        if c['bigg_id'] != None:
            for bid in c['bigg_id']:
                if Metabolite(bid + '_c') in model.metabolites:
                    # meta_dict[model.metabolites.get_by_id(bid+'_c')] = -1
                    meta_dict[model.metabolites.get_by_id(bid + '_c')] = rxn_coff[cpd]
                    break

        else:
            # print('Something is error. Mannully add compound %s to metabolites ' %(cpd))
            c_name = c['kegg_id'] + '_' + str(c['pubchem'])
            names[c_name] = Metabolite(id=cpd + '_c',
                                       compartment='c',
                                       name=c_name,
                                       formula=c['formula'])
            # meta_dict[names.get(c_name)] = -1
            meta_dict[names.get(c_name)] = rxn_coff[cpd]

    for cpd in products:
        cpd = re.search(r'.*(C\d{5}).*', cpd).group(1)
        # c = MyCompound(cpd)
        c = cpd_dict[cpd]

        if c['bigg_id'] != None:
            for bid in c['bigg_id']:
                if Metabolite(bid + '_c') in model.metabolites:
                    # meta_dict[model.metabolites.get_by_id(bid+'_c')] = 1
                    meta_dict[model.metabolites.get_by_id(bid + '_c')] = rxn_coff[cpd]
                    break

        else:
            # print('Something is error. Mannully add compound %s to metabolites ' %(cpd))
            c_name = c['kegg_id'] + '_' + str(c['pubchem'])
            names[c_name] = Metabolite(id=cpd + '_c',
                                       compartment='c',
                                       name=c_name,
                                       formula=c['formula'])
            # meta_dict[names.get(c_name)] = 1
            meta_dict[names.get(c_name)] = rxn_coff[cpd]

    return meta_dict


organism_list = ["ecoli", "yeast"]

def host(organism):
    if organism == organism_list[0]:
        dir = cfg['file_path'][organism+'_dir']
    if organism == organism_list[1]:
        dir = cfg['file_path'][organism+'_dir']
    return Path(dir)


def host_model(dir, model_params):
    model_dir = os.path.abspath(str(dir) + '/' + model_params)
    # Check if the file exists
    if os.path.isfile(model_dir):
        return model_dir
    # Return None if the file does not exist
    return print("ERROR!")


def get_yield(model_dir, rxn_list):

    original = read_sbml_model(model_dir)
    model = original.copy()

    medium = model.medium
    medium["EX_o2_e"] = 20.0
    medium["EX_glc__D_e"] = 20.0
    model.medium = medium

    wt_growth = model.optimize()
    max_growth = wt_growth.objective_value
    min_growth = 0.8 * max_growth

    rxn_names = locals()

    obj_rxn_id = None

    for rxn in rxn_list:

        # r = MyReaction(rxn)
        if rxn not in rxn_dict.keys():
            return 0

        r = rxn_dict[rxn]

        # First Judge rxn is already in model.reactions
        if r['bigg_id'] != None:
            for rid in r['bigg_id']:
                if Reaction(rid) in model.reactions:
                    # print('%s is already in model' %(rxn))
                    obj_rxn_id = rid
                    break

        else:
            rxn_names[rxn] = cobra.Reaction(rxn)
            rxn_names[rxn].lower_bound = -1000
            rxn_names[rxn].upper_bound = 1000

            rxn_names[rxn].add_metabolites(add_cpd2rxn(model, rxn))

            model.add_reactions([rxn_names[rxn]])

            # print('%s is done!' %(rxn))
            obj_rxn_id = rxn

    # model.add_boundary(model.metabolites.get_by_id(), type='demand')

    with model:
        medium = model.medium

        model.objective = model.reactions.get_by_id(obj_rxn_id)
        solution = model.optimize()
        max_biomass = solution.objective_value

        if max_biomass < min_growth:
            print("This pathway is not feasible!")
            return -999

        else:
            if model.reactions.get_by_id('EX_glc__D_e').flux == 0:
                print("This pathway cannot get EX_glc__D_e, Maximum theoretical yield = 0")
                return 0
            else:
                maximum_yield = max_biomass / (
                            -1 * (model.reactions.get_by_id('EX_glc__D_e').flux))  # Target production[mmol/gDW*h]
                # print('Maximum productivity =', max_biomass, 'mmol/gDW*h')

                # if maximum_yield > 49:
                #     print("This pathway EX_glc__D_e is 50, so the Maximum theoretical yield = 0")
                #    return 0.01
                #else:
                    #print('Maximum theoretical yield =', maximum_yield, 'mmol-Product/mmol-Glucose')
                return maximum_yield


# def main1():
#     # rxn_list = ['R07265', 'R11306', 'R04411', 'R04410', 'R09127', 'R06973']
#     rxn_list = ['R01788', 'R02737', 'R02780', 'R02628', 'R00724']
#
#     # org = host("ecoli")
#     # #model_dir = host_model(org, 'iJO1366.xml')  # 10.546250000000043
#     # model_dir = host_model(org, 'iEcHS_1320.xml')  # 10.546250000000013
#     org = host("yeast")
#     model_dir = host_model(org, 'iMM904.xml')  # 3.2750000000002366
#
#     print(get_yield(model_dir, rxn_list))

def maximum_yield(rxn_list, organism, model):

    org = host(organism)
    model_dir = host_model(org, model)  # 10.546250000000043
    #print(get_yield(model_dir, rxn_list))
    get_yield(model_dir, rxn_list)
    #print(yields)

if __name__ == "__main__":

    kegg_pathway = ['R01788', 'R02737', 'R02780', 'R02628', 'R00724']
    organism = "ecoli"
    model = "iJO1366.xml"
    maximum_yield(kegg_pathway, organism, model)
#
#     org = host(organism)
#     model_dir = host_model(org, model)
#     yields = get_yield(model_dir, kegg_pathway)
#     print(yields)



