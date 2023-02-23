"""query: Module for executing SPARQL queries via Rhea.

author(s): @christinehc
"""
# -------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------

from http.client import IncompleteRead
from typing import Optional

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper import sparql_dataframe as spdf

# -------------------------------------------------------------------
# Constants, URLs, and Parameters
# -------------------------------------------------------------------

UNIPROT = "https://sparql.uniprot.org/sparql/"
RHEA = "https://sparql.rhea-db.org/sparql"

# -------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------


def clean_query(query: str) -> str:
    """Format query cleanly for submission with valid newlines.

    Parameters
    ----------
    query : str
        Input query

    Returns
    -------
    str
        Formatted query

    """
    # compress whitespaces while preserving newlines
    return "\n".join(" ".join(line.split()) for line in query.splitlines())


def uniprot_query(identifier: str) -> str:
    """Write query for Rhea rxn IDs associated with UniProt ID.

    Parameters
    ----------
    identifier : str
        UniProt identifier

    Returns
    -------
    str
        Query for Rhea SPARQL endpoint

    """
    query = """
    PREFIX rh: <http://rdf.rhea-db.org/>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
    PREFIX chebi: <http://purl.obolibrary.org/obo/chebi/>

    SELECT DISTINCT ?uniprot ?identifier ?chebiId ?name ?smiles (group_concat(?rheaId; separator=',') as ?rheaIds)
    WHERE
    {{
        SERVICE <https://sparql.uniprot.org/sparql>
        {{
            ?uniprot up:annotation/up:catalyticActivity/up:catalyzedReaction ?rhea .
            VALUES (?uniprot) {{ (uniprotkb:{}) }}
        }}

        ?rhea rh:accession ?accession .
        ?rhea rh:equation ?equation .
        ?rhea rh:side/rh:contains/rh:compound ?compound .
        ?compound (rh:chebi|(rh:reactivePart/rh:chebi)|(rh:underlyingChebi/rh:chebi)) ?chebi .
        ?compound rh:name ?name .
        ?chebi chebi:smiles ?smiles .
        BIND(strafter(str(?uniprot),"http://purl.uniprot.org/uniprot/") as ?identifier)
        BIND(strafter(str(?chebi),"http://purl.obolibrary.org/obo/") as ?chebiId)
        BIND(strafter(str(?accession),"RHEA:") as ?rheaId)
    }}
    GROUP BY ?uniprot ?identifier ?chebiId ?name ?smiles
    """.format(
        identifier
    )
    return clean_query(query)


def ec_query(identifier: str) -> str:
    """Write query for Rhea rxn IDs associated with EC ID.

    Parameters
    ----------
    identifier : str
        EC identifier

    Returns
    -------
    str
        Query for Rhea SPARQL endpoint

    """
    identifier = f"^{identifier}$"
    query = """
    PREFIX rh: <http://rdf.rhea-db.org/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX ec: <http://purl.uniprot.org/enzyme/>
    PREFIX chebi: <http://purl.obolibrary.org/obo/chebi/>

    SELECT ?ec ?identifier ?chebiId ?name ?smiles (group_concat(?rheaId; separator=',') as ?rheaIds)
    WHERE
    {{
      ?rhea rdfs:subClassOf rh:Reaction .
      ?rhea rh:accession ?accession .
      ?rhea rh:ec ?ec .

      BIND(strafter(str(?ec),str(ec:)) as ?identifier)
      FILTER (regex(?identifier,"{}")) 

      ?rhea rh:equation ?equation .
      ?rhea rh:side/rh:contains/rh:compound ?compound .
      ?compound (rh:chebi|(rh:reactivePart/rh:chebi)|(rh:underlyingChebi/rh:chebi)) ?chebi .
      ?compound rh:name ?name .
      ?chebi chebi:smiles ?smiles .
      BIND(strafter(str(?chebi),"http://purl.obolibrary.org/obo/") as ?chebiId)
      BIND(strafter(str(?accession),"RHEA:") as ?rheaId)
    }}
    GROUP BY ?ec ?identifier ?chebiId ?name ?smiles
    """.format(
        identifier
    )
    return clean_query(query)


def rhea_query(identifier: str) -> str:
    """Write query for SMILES associated with RHEA ID.

    Parameters
    ----------
    identifier : str
        RHEA identifier

    Returns
    -------
    str
        Query for Rhea SPARQL endpoint

    """
    query = """
    PREFIX rh: <http://rdf.rhea-db.org/>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
    PREFIX chebi: <http://purl.obolibrary.org/obo/chebi/>

    SELECT DISTINCT ?chebiId ?name ?smiles (group_concat(?rheaId; separator=',') as ?rheaIds)
    WHERE
    {{
      VALUES (?rhea){ { (rh:{}) }}
      ?rhea rh:accession ?accession .
      ?rhea rh:equation ?equation .
      ?rhea rh:side/rh:contains/rh:compound ?compound .
      ?compound (rh:chebi|(rh:reactivePart/rh:chebi)|(rh:underlyingChebi/rh:chebi)) ?chebi .
      ?compound rh:name ?name .
      ?chebi chebi:smiles ?smiles .
      BIND(strafter(str(?chebi),"http://purl.obolibrary.org/obo/") as ?chebiId)
      BIND(strafter(str(?accession),"RHEA:") as ?rheaId)
    }}
    GROUP BY ?chebiId ?name ?smiles
    """.format(
        identifier
    )
    return clean_query(query)


def id2smiles(
    identifier: str, endpoint: str = RHEA, id_type: str = "uniprot"
) -> pd.DataFrame:
    """Convert database ID to SMILES via SPARQL query.

    Parameters
    ----------
    identifier : str
        Value of identifier
    endpoint : str, optional
        SPARQL endpoint, one of ["UNIPROT", "RHEA"]
        by default UNIPROT
    id_type : str, optional
        Identifier type, by default "uniprot"
        Accepts "uniprot" / "up" for UniProt ID
                "enzyme" / "ec" for EC ID
                "rhea" / "rh" for RHEA reaction ID

    Returns
    -------
    pd.DataFrame
        Dataframe with ChEBI, and SMILES matching input identifier

    Raises
    ------
    ValueError
        Raised for invalid `id_type` options.
    """
    id_type = id_type.lower()
    if id_type in ["uniprot", "up"]:
        query = uniprot_query(identifier).encode("utf-8")
    elif id_type in ["enzyme", "ec"]:
        query = ec_query(identifier).encode("utf-8")
    elif id_type in ["rhea", "rh"]:
        query = rhea_query(identifier).encode("utf-8")
    else:
        raise ValueError("Parameter `id_type` must be 'uniprot', 'ec', or 'rhea'.")
    df = spdf.get_sparql_dataframe(endpoint, query)
    df["chebiId"] = [c.replace("CHEBI_", "CHEBI:") for c in df["chebiId"]]
    return df


def chebi2property(chebi: str, prop: str = "SMILES") -> Optional[str]:
    """Search ChEBI for property by ChEBI ID.

    Parameters
    ----------
    chebi : str
        ChEBI identifier
    prop : str, optional
        Name of property, by default "SMILES"

    Returns
    -------
    Optional[str]
        Value of queried property; returns None if none found
    """
    if "CHEBI" not in chebi:
        chebi = f"CHEBI:{chebi}"
    url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi}"
    tabs = list()
    try:
        tabs = [t for t in pd.read_html(url) if prop in t.values]
    except IncompleteRead:
        pass
    if len(tabs) > 0:
        return tabs[0].values[0][1]
    return


def chebis2smilesdf(chebis: list[str]) -> pd.DataFrame:
    """Populate dataframe with compound info given ChEBI IDs.

    Parameters
    ----------
    chebis : list[str]
        List of ChEBI IDs

    Returns
    -------
    pd.DataFrame
        Dataframe containing ChEBIs, names, and SMILES
    """
    smiles = [chebi2property(c, "SMILES") for c in chebis]
    names = [chebi2property(c, "ChEBI Name") for c in chebis]
    return pd.DataFrame({"chebiId": chebis, "name": names, "smiles": smiles})
