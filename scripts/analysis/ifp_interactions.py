from mol.entry import *
import napoli
#from napoli import *
from util.default_values import *

#from mol.entry import recover_entries_from_entity
from MyBio.util import download_pdb
from MyBio.PDB.PDBParser import PDBParser
#from mol.entry import *
from mol.interaction.view import InteractionViewer
from mol.interaction.fp.view import ShellViewer

from mol.interaction.filter import InteractionFilter
from mol.interaction.calc import InteractionCalculator
import rdkit.Chem
import sys,os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import numpy as np
#import rdkit.Chem.Descriptors as descr

def RetrieveMol2Block(fileLikeObject, delimiter="@<TRIPOS>MOLECULE"):
    """generator which retrieves one mol2 block at a time
    """
    mol2 = []
    for line in fileLikeObject:
        if line.startswith(delimiter) and mol2:
            yield "".join(mol2)
            mol2 = []
        mol2.append(line)
    if mol2:
        yield "".join(mol2)

def count_unsatisfied_donors(entry):
    # check unsatisfied hbond donor
    count_unsatisfied = 0
    for atm in entry['donor_atm_list']:
        flag_satisfied = False
        if not ('Hydrogen bond' in entry['interactions'].keys()):
            #flag_filters = False
            print("%s doesn't have hbond at the binding site" % (entry['id']))
            break
        for hbond_type in entry['interactions']['Hydrogen bond']:
            #print(atm)
            #print(hbond_type)
            if ((atm[0] in hbond_type[0]) or (atm[0] in hbond_type[1])):
                flag_satisfied = True
                #print(atm)
                #print(hbond_type)
                break
        if flag_satisfied == False:
            count_unsatisfied = count_unsatisfied + 1
    return count_unsatisfied

def count_unsatisfied_acceptors(entry):
    # check unsatisfied hbond acceptor
    count_unsatisfied = 0
    for atm in entry['acceptor_atm_list']:
        flag_satisfied = False
        if not ('Hydrogen bond' in entry['interactions'].keys()):
            #flag_filters = False
            print("%s doesn't have hbond at the binding site" % (entry['id']))
            break
        for hbond_type in entry['interactions']['Hydrogen bond']:
            #print(atm)
            #print(hbond_type)
            if ((atm[0] in hbond_type[0]) or (atm[0] in hbond_type[1])):
                flag_satisfied = True
                #print(atm)
                #print(hbond_type)
                break
        if flag_satisfied == False:
            count_unsatisfied = count_unsatisfied + 1
    return count_unsatisfied

def find_inter(entry, filter):
    flag_filter = False
    # break the loop if the entry doesn't have this interaction type
    try:
        if filter[0] == 'Ionic':
            for type in entry['interactions'][filter[0]]:
                if (filter[1] in type[0]):
                    # check if they are in hydrogen bound distances
                    for hbond_type in entry['interactions']['Hydrogen bond']:
                        if ((filter[1] in hbond_type[0]) and (type[1] in hbond_type[1])) or ((filter[1] in hbond_type[1]) and (type[1] in hbond_type[0])):
                            flag_filter = True
                elif (filter[1] in type[1]):
                    # check if they are in hydrogen bound distances
                    for hbond_type in entry['interactions']['Hydrogen bond']:
                        if ((filter[1] in hbond_type[0]) and (type[0] in hbond_type[1])) or ((filter[1] in hbond_type[1]) and (type[0] in hbond_type[0])):
                            flag_filter = True
        else: 
            for type in entry['interactions'][filter[0]]:
                if (filter[1] in type[0]) or (filter[1] in type[1]):
                    # interaction found
                    flag_filter = True
    except:
        print("%s %s doesn't exist in the docked pose of %s\n" % (filter[0],filter[1],entry['id']))
    return flag_filter

if len(sys.argv) != 6:
    print("%s <sdf_file> <PDB_ID> <PDB_path> <fp_len> <output_profix>" % sys.argv[0])
    exit()

sdf_filename  = sys.argv[1]
PDB_ID        = sys.argv[2]
PDB_path      = sys.argv[3]
fp_len        = int(sys.argv[4])
output_profix = sys.argv[5]
#working_path  = sys.argv[5]


mols = []
#suppl = rdkit.Chem.SDMolSupplier(sdf_filename)
suppl = rdkit.Chem.SDMolSupplier(sdf_filename,removeHs=False)
for mol in suppl:
    try:
        name_ori = mol.GetProp('_Name')
        name = name_ori.split()[0]
        mol.SetProp('_Name',name)
        #print(mol.GetProp('_Name'))
        mols.append(mol)
    except:
        continue
count = 0
entries = []
for mol in mols:
    entry = MolEntry(pdb_id=PDB_ID, mol_id=mol.GetProp('_Name'), mol_obj=mol, mol_obj_type='rdkit', mol_file="./poses_extact_for_getposes_parallel_"+output_profix+".sdf")
    #print(entry)
    entries.append(entry)

ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter(),strict_donor_rules=True)
#ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter())

#working_path = "%s/tmp/local_proj" % NAPOLI_PATH

print("Number of entries to be processed: %d." % len(entries))

opt = {}
opt["working_path"] = PDB_path
opt["overwrite_path"] = False
opt["entries"] = entries
opt["preload_mol_files"] = True
opt["inter_calc"] = ic
opt["mol_obj_type"] = 'rdkit'
#opt["try_h_addition"] = True
# if you have hydrogen-added pdb and sdf files, don't turn it on!!!
opt["try_h_addition"] = False
opt["amend_mol"] = True
#opt["amend_mol"] = False
#from napoli import PLI_PROJECT
#project_type = PLI_PROJECT
#project_type = napoli.PLI_PROJECT
#opt["project_type"] = project_type
opt["ifp_output"] = "poses_"+output_profix+".csv"
#opt["add_non_cov"] = True
#opt["add_atom_atom"] = True
#opt["add_proximal"] = False
opt["ifp_num_levels"] = 5
opt["ifp_radius_step"] = 1
opt["ifp_length"] = fp_len
#opt["mol_obj_type"] = 'rdkit'
opt["calc_mfp"] = False
opt["calc_ifp"] = True
opt["pdb_path"] = PDB_path

from napoli import InteractionProject
proj_obj = InteractionProject(**opt)
proj_obj()
results = proj_obj.result
#filters = [['Ionic','ASP-115']]
#filters = [['Ionic','ASP-115'],['Hydrogen bond','ASP-115']]
#filters = [['Ionic','ASP-155'],['Hydrogen bond','THR-160'],['Hydrogen bond','SER-239'],['Hydrogen bond','SER-242'],['Van der Waals','GLY-369']]
# HB to HIS163, GLY143, GLU166 and hydrophobic pocket
#filters = [['Hydrogen bond','HIS-163'],['Hydrogen bond','GLY-143'],['Hydrogen bond','GLU-166'],['Hydrophobic','ASP-187'],['Van der Waals','MET-49'],['Hydrophobic','MET-49']]
#filters = [['Hydrogen bond','ARG-88'],['Hydrogen bond','ARG-322']]
filters = [['Hydrogen bond','GLY-333'],['Hydrogen bond','ALA-353'],['Hydrogen bond','TYR-368']]
outputfile = open(output_profix+".interactions.csv",'w')

# write the header
# id #_donors #_acceptors #_unstatisfied_donors #_unstatisfied_acceptors interaction_residue_id
outputfile.write("id,#_donors,#_acceptors,#_unstatisfied_donors,#_unstatisfied_acceptors")
for filter in filters:
    outputfile.write(",%s/%s" % (filter[0],filter[1]))
outputfile.write("\n")

for entry in results:
    # zincid
    id = entry['id'].split(':')[1]
    #print(id)
    #print(entry['interactions']['Hydrogen bond'])
    #print(entry['interactions'])
    # #_donors
    #print(entry['donor_atm_list'])
    num_donors = len(entry['donor_atm_list'])
    # #_acceptors
    #print(entry['acceptor_atm_list'])
    num_acceptors = len(entry['acceptor_atm_list'])
    # #_unstatisfied_donors
    num_unstatisfied_donors = count_unsatisfied_donors(entry)
    # #_unstatisfied_acceptors
    num_unstatisfied_acceptors = count_unsatisfied_acceptors(entry)
    # interactions
    interaction_flags = []
    num_interactions = len(filters)
    for filter in filters:
        flag = find_inter(entry, filter)
        if flag == True:
            interaction_flags.append(1)
        else:
            interaction_flags.append(0)
    if not (num_interactions == len(interaction_flags)):
        print("num_interactions(%d) != interaction_flags(%d)" % (num_interactions, len(interaction_flags)))
        exit()
    outputfile.write("%s,%d,%d,%d,%d" % (id,num_donors,num_acceptors,num_unstatisfied_donors,num_unstatisfied_acceptors))
    for flag in interaction_flags:
        outputfile.write(",%d" % flag)
    outputfile.write("\n")

print("DONE!!!")
