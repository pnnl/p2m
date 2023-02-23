"""utils: Define general utility functions.

authors: @brykpnl
"""

import re
import libchebipy
import numpy as np
from numpy.typing import NDArray
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


def chebi_rels(chebi_id: str, wait: float = 0.1) -> list[str]:
    """Send external requests for related ChEBI molecules.

    Parameters
    ----------
    chebi_id : str
        _description_
    wait : float, optional
        _description_, by default 0.1

    Returns
    -------
    _type_
        _description_
    """
    chebi = libchebipy.ChebiEntity(chebi_id)
    chebi_rels = chebi.get_incomings()
    chebis = [z._Relation__target_chebi_id for z in chebi_rels]
    time.sleep(wait)
    return chebis


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


def query_external_uniprot(
    uids: list[str], n_reqs: int = 3, ext: str = "txt"
) -> NDArray:
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
    dfc["cleanedSmiles"] = dfc["smiles"].apply(
        standardize_smiles, args=(TautomerEnumerator(), Uncharger())
    )
    return dfc
