"""annotation: Annotate chemical identifiers with structural information.

authors: @brykpnl, @christinehc
"""
import os

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
        self, ids: list, ids_type: str = "uniprot", external_uniprot: bool = False
    ):
        """Run SPARQL query for specified IDs."""
        self.ids = ids
        self.ids_type = ids_type

        result = list()
        for id_ in ids:
            tmp = query.id2smiles(id_, id_type=ids_type)
            result.append(tmp)
        result = pd.concat(result, ignore_index=True)

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
