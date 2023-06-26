"""annotation: Annotate chemical identifiers with structural information.

authors: @brykpnl, @christinehc
"""
import os
from time import sleep

import pandas as pd
from . import query, utils


class Annotation:
    """Class for mapping UniProt or EC identifiers to metabolites."""

    def __init__(self):
        self.ids = []
        self.ids_type = None
        self.rhea_ids = []
        self.rxns = []

    def query_ids(
        self,
        ids: list,
        ids_type: str = "uniprot",
        external_uniprot: bool = False,
        sleep_time: float = 3.0,
    ):
        """Run SPARQL query for specified IDs.

        Parameters
        ----------
        ids : list
            List of UniProt/EC identifiers
        ids_type : str, optional
            Identifier type, by default "uniprot"
            Accepts "uniprot" / "up" for UniProt ID
                    "enzyme" / "ec" for EC ID
                    "rhea" / "rh" for RHEA reaction ID
        external_uniprot : bool, optional
            If True, applies UniProt REST for IDs missing rxn info,
                by default False
        sleep_time : float, optional
            ait time (s) between query requests, by default 3.0
        """
        self.ids = ids
        self.ids_type = ids_type

        result = list()
        for id_ in ids:
            tmp = query.id2smiles(id_, id_type=ids_type)
            sleep(sleep_time)
            result.append(tmp)
        result = pd.concat(result, ignore_index=True)

        if result.empty:
            raise ValueError(
                "Queried IDs return no results. Check that input"
                " identifiers are annotated in the Rhea database"
                " and the correct `--type` flag was indicated."
            )

        self.rhea_ids.extend(
            list(set([i for item in result["rheaIds"].values for i in item.split(",")]))
        )
        self.missing_uni_ids = [
            i for i in ids if i not in result["identifier"].unique()
        ]

        if external_uniprot:
            self.external_uniprot_rhea = utils.query_external_uniprot(
                self.missing_uni_ids
            )
            additional_rxns = [
                n for n in self.external_uniprot_rhea if n not in self.rhea_ids
            ]
            self.rhea_ids.extend(additional_rxns)  # append additional IDs

            additional_df = list()
            if len(additional_rxns) > 0:
                for id_ in additional_rxns:
                    additional_df.append(query.id2smiles(id_, id_type="rhea"))
            additional_df = pd.concat(additional_df, ignore_index=True)
            result = pd.concat([result, additional_df], ignore_index=True)

        self.smiles = result.drop_duplicates(subset="smiles", ignore_index=True)
        self.smiles_stars = self.smiles[self.smiles.smiles.str.contains("\*")]
        self.smiles_nostars = self.smiles[~self.smiles.smiles.str.contains("\*")]

    def star_to_smiles(self):
        """Externally query ChEBI database for related compounds"""
        related_chebis = [
            utils.chebi_rels(chebi)
            for chebi in self.smiles_stars["chebiId"]
            if "CHEBI" in chebi
        ]
        flattened_related = list(
            set(["CHEBI:" + str(i) for s in related_chebis for i in s])
        )
        matched_chebis = query.chebis2smilesdf(flattened_related)
        all_smiles = matched_chebis[matched_chebis["smiles"].notna()]
        nostars = pd.concat(
            [self.smiles_nostars, all_smiles[~all_smiles.smiles.str.contains("\*")]]
        )
        self.all_nostars = nostars.drop_duplicates(subset="chebiId", ignore_index=True)

    def clean_smiles(self):
        """
        Pass all SMILES strings in the dataframe through standardization steps
        that desalt, uncharge, and generate a canonical tautomer for each
        SMILES entry in the dataframe.
        """
        self.clean_smiles_nostars = utils.clean_smiles_df(self.smiles_nostars)
        if hasattr(self, "all_nostars"):
            self.clean_all_nostars = utils.clean_smiles_df(self.all_nostars)

    def export_smiles(self, export_path: str):
        """Export SMILES dataframes with ChEBI identifiers to appropriate locationsin the output folder."""
        self.smiles_stars.to_csv(
            os.path.join(export_path, "smiles_stars.tsv"), index=False, sep="\t"
        )
        if hasattr(self, "clean_smiles_nostars"):
            self.clean_smiles_nostars.to_csv(
                os.path.join(export_path, "clean_smiles_nostars.tsv"),
                index=False,
                sep="\t",
            )
            if hasattr(self, "clean_all_nostars"):
                self.clean_all_nostars.to_csv(
                    os.path.join(export_path, "clean_all_nostars.tsv"),
                    index=False,
                    sep="\t",
                )
        else:
            self.smiles_nostars.to_csv(
                os.path.join(export_path, "smiles_nostars.tsv"), index=False, sep="\t"
            )
            if hasattr(self, "all_nostars"):
                self.all_nostars.to_csv(
                    os.path.join(export_path, "smiles_all_nostars.tsv"),
                    index=False,
                    sep="\t",
                )
