import re
import libchebipy
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.MolStandardize.rdMolStandardize import (
    TautomerEnumerator,
    Uncharger,
    StandardizeSmiles,
)
from openbabel import pybel
import urllib.request
import time
import ssl
import logging
import pandas as pd

UNIPROT_REST_URL = "https://rest.uniprot.org/uniprotkb/"


def get_master_rhea(rhea_map, ids):
    """Identify which rhea ID's are found in the local database and those that are not."""
    rhea_ids = list(rhea_map["ID"])
    matched_ids = list(set(rhea_ids).intersection(set(ids)))
    missing_ids = list(set(ids) - set(rhea_ids))
    found_entries = rhea_map.loc[rhea_map["ID"].isin(matched_ids)]
    found_rheas = list(found_entries["MASTER_ID"])
    return found_rheas, missing_ids


def get_rxn_ids(rhea_ids, rhea_directions):
    """Read in ambiguous Rhea identifiers.
    Return a dict containing rhea identifiers with local IDs
    with those not in rhea_directions.
    """
    master_rxns = rhea_directions.loc[rhea_directions["RHEA_ID_MASTER"].isin(rhea_ids)]
    master_found = list(master_rxns["RHEA_ID_MASTER"])
    lr_master = list(master_rxns["RHEA_ID_LR"])
    rl_master = list(master_rxns["RHEA_ID_RL"])

    lr_only = list(
        rhea_directions.loc[rhea_directions["RHEA_ID_LR"].isin(rhea_ids)]["RHEA_ID_LR"]
    )
    rl_only = list(
        rhea_directions.loc[rhea_directions["RHEA_ID_RL"].isin(rhea_ids)]["RHEA_ID_RL"]
    )
    bi_lr = list(
        rhea_directions.loc[rhea_directions["RHEA_ID_BI"].isin(rhea_ids)]["RHEA_ID_LR"]
    )
    bi_rl = list(
        rhea_directions.loc[rhea_directions["RHEA_ID_BI"].isin(rhea_ids)]["RHEA_ID_RL"]
    )

    bi_both = bi_lr + bi_rl
    all_found_rhea = lr_master + rl_master + lr_only + rl_only + bi_both
    not_found_rhea = list(
        set(rhea_ids)
        - set(master_found)
        - set(lr_master)
        - set(rl_master)
        - set(bi_both)
    )
    rhea_rxn_ids = {"local": all_found_rhea, "external": not_found_rhea}
    return rhea_rxn_ids


def chebi_rels(chebi_id, wait=0.1):
    """Send external requests for related ChEBI molecules."""
    chebi = libchebipy.ChebiEntity(chebi_id)
    chebi_rels = chebi.get_incomings()
    chebis = [z._Relation__target_chebi_id for z in chebi_rels]
    time.sleep(wait)
    return chebis


def get_chebi_smiles(chebi_frame: pd.DataFrame, chebis: list[str]) -> pd.DataFrame:
    """Get SMILES strings from ChEBI data frame."""
    entries = chebi_frame.loc[chebi_frame["ChEBI_ID"].isin(chebis)]
    entries.reset_index(inplace=True, drop=True)
    return entries


def external_rhea(rhea_ids: list[str], wait: float = 0.1) -> list:
    """Send external requests to RHEA db.

    Parameters
    ----------
    rhea_ids : list[str]
        List of RHEA reaction identifiers to be queried
    wait : float, optional
        Wait time between requests in seconds, by default 0.1

    Returns
    -------
    list
        List of
    """
    rxns = list()
    for rhea in rhea_ids:
        try:
            query = "https://www.rhea-db.org/rest/1.0/ws/reaction/rxn/" + rhea
            contents = urllib.request.urlopen(query, context=ssl.SSLContext()).read()
            rxns.append(rdChemReactions.ReactionFromRxnBlock(contents))
        except Exception:
            logging.info(
                "Could not retrieve guessed external rhea reaction: {}".format(rhea)
            )
        time.sleep(wait)
    return rxns


def convert_rxn(rxn_path):
    """Convert .rxn file to an RDKit reaction."""
    try:
        return rdChemReactions.ReactionFromRxnFile(rxn_path)
    except Exception:
        return "ERROR: {}".format(rxn_path)


def rxn_to_smiles(rxn):
    """Obtain SMILES strings of produts from RDKit reactions."""
    try:
        products = rxn.GetProducts()
        all_smiles = {mol.GetProp("_Name"): Chem.MolToSmiles(mol) for mol in products}
        return all_smiles
    except Exception:
        return "ERROR: {}".format(rxn)


def get_chebi_names(chebi_path, smiles_df):
    """Get ChEBI molecule names from ids."""
    chebi_names = pd.read_csv(chebi_path, sep="\t", header=None)
    chebi_names.columns = ["ChEBI_ID", "Name"]
    found_names = chebi_names.loc[chebi_names["ChEBI_ID"].isin(smiles_df["ChEBI_ID"])]
    full_df = pd.merge(found_names, smiles_df, on="ChEBI_ID", how="outer")
    return full_df


def query_external_uniprot(
    uids: list[str], n_reqs: int = 3, ext: str = "txt"
) -> dict[str, list]:
    """Query the UniProt database for Rhea reaction ids via REST.

    Parameters
    ----------
    uids : list[str]
        List of UniProt IDs
    n_reqs : int, optional
        Maximum requests to server, by default 3
    ext : str, optional
        Page extension for REST, by default "txt"
        (Options: "txt", "rdf")

    Returns
    -------
    dict[str, list]
        Mapping of UniProt ID to list of Rhea reaction IDs
    """
    url = UNIPROT_REST_URL

    mapping = {}
    for uid in uids:
        i = 0

        req = urllib.request.Request(url + f"{uid}.{ext}")
        while i < n_reqs:
            try:
                with urllib.request.urlopen(req) as f:
                    response = f.read()
                break
            except Exception:
                logging.info("Failed UniProt request {}. Repeating...".format(i))
                i += 1
        rd = response.decode("utf-8")
        selector = r"RHEA:([0-9]*)"
        rhea_selector = re.compile(selector)
        rxn_ids = re.findall(rhea_selector, rd)
        mapping[uid] = rxn_ids
    return np.unique(np.hstack(list(mapping.values())))


def standardize_smiles(smiles, te, uc):
    """Standardize a SMILES string by desalting, uncharging, and canonicalizing."""
    try:
        mol = pybel.readstring("smi", smiles)
        mol.OBMol.StripSalts()
        desalt_mol = Chem.MolFromSmiles(mol.write("can").strip())
    except Exception:
        logging.info(
            "Couldn't desalt SMILES:\n{}\nWill still try to standardize...\n".format(
                smiles
            )
        )
        try:
            desalt_mol = Chem.MolFromSmiles(smiles)
        except Exception:
            logging.info("Couldn't desalt {}. Skipping.".format(smiles))
            return None
    try:
        standardized = Chem.MolToSmiles(te.Canonicalize(uc.uncharge(desalt_mol)))
        return standardized
    except Exception:
        logging.info("Couldn't standardize {}. Skipping.".format(smiles))
        return None


def clean_smiles_df(df):
    """Standardize SMILES strings."""
    handle = pybel.ob.OBMessageHandler()
    handle.SetOutputLevel(0)
    dfc = df.copy()
    dfc["cleaned_smiles"] = dfc["smiles"].apply(
        standardize_smiles, args=(TautomerEnumerator(), Uncharger())
    )
    return dfc


def get_chebi_info(mol):
    """Get ChEBI id, name, and smiles strings from mol file."""
    chebi_id = mol.GetProp("ChEBI ID")
    chebi_name = mol.GetProp("ChEBI Name")
    smiles = Chem.MolToSmiles(mol)
    return [chebi_id, chebi_name, smiles]
