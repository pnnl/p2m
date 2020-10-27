import re
import libchebipy
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator, Uncharger, StandardizeSmiles
import pybel
import urllib.request
import time
import ssl
import pandas as pd

def get_master_rhea(rhea_map, ids):
    rhea_ids = list(rhea_map['ID'])
    matched_ids = list(set(rhea_ids).intersection(set(ids))) # All Uniprot identifiers that were found with a local cross-reference
    missing_ids = list(set(ids)-set(rhea_ids)) # All Uniprot identifiers not found with local cross-reference
    found_entries = rhea_map.loc[rhea_map['ID'].isin(matched_ids)]
    found_rheas = list(found_entries['MASTER_ID'])
    return(found_rheas, missing_ids)

#TODO: Return only those based on argument passed. Right now, it's everything.
#As of now, everything is selected (both directions). 
#TODO: Also clean up as much as possible

def get_rxn_ids(rhea_ids, rhea_directions):
    """
    Read in ambiguous Rhea identifiers and return a dict containing rhea identifiers with local ID's along
    with those that aren't in rhea_directions.
    """
    master_rxns = rhea_directions.loc[rhea_directions['RHEA_ID_MASTER'].isin(rhea_ids)]
    master_found = list(master_rxns['RHEA_ID_MASTER'])
    lr_master = list(master_rxns['RHEA_ID_LR'])
    rl_master = list(master_rxns['RHEA_ID_RL'])
    lr_only = list(rhea_directions.loc[rhea_directions['RHEA_ID_LR'].isin(rhea_ids)]['RHEA_ID_LR'])
    rl_only = list(rhea_directions.loc[rhea_directions['RHEA_ID_RL'].isin(rhea_ids)]['RHEA_ID_RL'])
    bi_lr = list(rhea_directions.loc[rhea_directions['RHEA_ID_BI'].isin(rhea_ids)]['RHEA_ID_LR'])
    bi_rl = list(rhea_directions.loc[rhea_directions['RHEA_ID_BI'].isin(rhea_ids)]['RHEA_ID_RL'])
    bi_both = bi_lr + bi_rl
    all_found_rhea = lr_master + rl_master + lr_only + rl_only + bi_both
    not_found_rhea = list(set(rhea_ids) - set(master_found) - set(lr_master) - set(rl_master) - set(bi_both))
    rhea_rxn_ids = {'local': all_found_rhea, 'external': not_found_rhea}
    return(rhea_rxn_ids)

def chebi_rels(chebi_id, wait=.1):
    chebi = libchebipy.ChebiEntity(chebi_id)
    chebi_rels = chebi.get_incomings()
    chebis = [z._Relation__target_chebi_id for z in chebi_rels]
    time.sleep(wait)
    return(chebis)

def get_chebi_smiles(chebi_frame, chebis):
    entries = chebi_frame.loc[chebi_frame['ChEBI_ID'].isin(chebis)]
    entries.reset_index(inplace=True, drop=True)
    return(entries)

def external_rhea(rhea_ids, wait=.1):
    rxns = list()
    for rhea in rhea_ids:
        try:
            query = "https://www.rhea-db.org/rest/1.0/ws/reaction/rxn/" + rhea
            contents = urllib.request.urlopen(query, context=ssl.SSLContext()).read()
            rxns.append(rdChemReactions.ReactionFromRxnBlock(contents))
        except Exception:
            print("Could not retrieve guessed external rhea reaction: {}".format(rhea))
        time.sleep(wait)
    return(rxns)
                
def convert_rxn(rxn_path):
    try:
        return rdChemReactions.ReactionFromRxnFile(rxn_path)
    except Exception:
        return "ERROR: {}".format(rxn_path)

def rxn_to_smiles(rxn):
    try:
        products = rxn.GetProducts()
        all_smiles = {mol.GetProp('_Name'): Chem.MolToSmiles(mol) for mol in products}
        return(all_smiles)
    except Exception:
        return("ERROR: {}".format(rxn))

def get_chebi_names(chebi_path, smiles_df):
    chebi_names = pd.read_csv(chebi_path, sep='\t', header=None)
    chebi_names.columns = ['ChEBI_ID', 'Name']
    found_names = chebi_names.loc[chebi_names['ChEBI_ID'].isin(smiles_df['ChEBI_ID'])]
    full_df = pd.merge(found_names, smiles_df, on='ChEBI_ID', how='outer')
    return(full_df)

def query_external_uniprot(uids, num_reqs=3):
    queries = " ".join(uids)
    params = {
    'from': 'ACC+ID',
    'to': 'ACC',
    'format': 'tab',
    'columns': 'rhea-id',
    'query': queries
    }
    url = 'https://www.uniprot.org/uploadlists/'
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    req.add_header('Content_Type', 'form-data') #Correct header for uniprot
    i=0
    while (i < num_reqs):
        try:
            with urllib.request.urlopen(req) as f:
               response = f.read()
            break
        except Exception:
            print("Failed UniProt request {}. Repeating...".format(i))
            i += 1
    rd = response.decode('utf-8')
    selector=r'RHEA:([0-9]*)'
    rhea_selector = re.compile(selector)
    rxn_ids = re.findall(rhea_selector, rd)
    return(rxn_ids)
    
def standardize_smiles(smiles, te, uc):
    """
    Standardize a SMILES string.
    Try-except blocks will cause the function to return None instead of halting if an error is encountered.
    """
    try:
        mol = pybel.readstring('smi', smiles)
        mol.OBMol.StripSalts()
        desalt_mol = Chem.MolFromSmiles(mol.write('can').strip())
    except Exception:
        print("Couldn't desalt SMILES:\n{}\nWill still try to standardize...\n".format(smiles))
        try:
            desalt_mol = Chem.MolFromSmiles(smiles)
        except Exception:
            return(None)
    try:
        standardized = Chem.MolToSmiles(te.Canonicalize(uc.uncharge(desalt_mol)))
        return(standardized)
    except Exception:
        return(None)
    
def clean_smiles_df(df):
    handle = pybel.ob.OBMessageHandler()
    handle.SetOutputLevel(0)
    dfc = df.copy()
    dfc['Cleaned_SMILES'] = dfc['SMILES'].apply(standardize_smiles, args=(TautomerEnumerator(), Uncharger()))
    return(dfc)
    
def get_chebi_info(mol):
    chebi_id = mol.GetProp('ChEBI ID')
    chebi_name = mol.GetProp('ChEBI Name')
    smiles = Chem.MolToSmiles(mol)
    return([chebi_id, chebi_name, smiles])



