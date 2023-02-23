"""main: Run main P2M program.

authors: @brykpnl, @christinehc
"""
import argparse
import logging
import os
import sys
import time

from rdkit import RDLogger

from . import annotation

RDLogger.DisableLog("rdApp.*")


class _HelpAction(argparse._HelpAction):
    def __call__(self, parser, namespace, values, option_string=None):
        # print main help
        print(parser.format_help())

        # retrieve subparsers from parser
        subparsers_actions = [
            action
            for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        ]
        # there will probably only be one subparser_action,
        # but better safe than sorry
        for subparsers_action in subparsers_actions:
            # get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                print(f"p2m {choice}")
                print(subparser.format_help())
        parser.exit()


def run(
    ids_path: str,
    id_type: str,
    output_path: str,
    get_related_chebis: bool = False,
    clean_smiles: bool = False,
):
    """Run p2m to map identifiers to SMILES strings.

    Parameters
    ----------
    ids_path : str
        Path to plaintext file containing UniProt/EC identifiers
    id_type : str
        Identifier type; select from ["uniprot", "up", "enzyme", "ec"]
    output_path : str
        Path where to save output files
    get_related_chebis : bool, optional
        _description_, by default False
    clean_smiles : bool, optional
        If True, cleans SMILES via RDKit, by default False
    """
    with open(ids_path) as f:
        ids = f.read().splitlines()

    annot = annotation.Annotation()

    logging.info("Mapping protein identifiers to Rhea identifiers...")
    annot.query_ids(ids, id_type)

    if get_related_chebis:
        logging.info(
            "Querying ChEBI externally for unspecified R-group containing compounds..."
        )
        annot.star_to_smiles()

    if clean_smiles:
        logging.info("Standardizing SMILES strings...")
        annot.clean_smiles()

    logging.info("Exporting output to: {}".format(output_path))
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    annot.export_smiles(output_path)


def main():
    """Run p2m to download required files or map identifiers to SMILES strings."""

    # main parser
    parser = argparse.ArgumentParser(
        description=(
            "---\n"
            "P2M\n"
            "---\n\n"
            "Identify metabolites (substrates and products) of proteins "
            "from UniProt identifiers or EC numbers. "
            "Requires internet access for external database queries. "
            "\n\nTo set up local databases for your first run,"
            " be sure to run the 'download' command."
        ),
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--help", "-h", action=_HelpAction, help="Show help documentation."
    )  # add custom help

    parser.add_argument(
        "--ids_path",
        "--input",
        "-i",
        type=str,
        help="Path to UniProt or EC identifiers; file must be a list of identifiers separated by newlines.",
        dest="ids_path",
        required=True,
    )
    parser.add_argument(
        "--ids_type",
        "--type",
        "-t",
        type=str.lower,
        help="Type of identifiers. One of: UniProt, EC",
        dest="ids_type",
        required=True,
    )
    parser.add_argument(
        "--output_path",
        "--output",
        "-o",
        type=str,
        help="Path to desired output folder.",
        dest="output_path",
        required=True,
    )
    parser.add_argument(
        "--complete_rgroups",
        "-r",
        action="store_true",
        help="Externally query the ChEBI database for substructure searches of compounds with R-groups. Default False.",
        dest="complete_rgroups",
    )
    parser.add_argument(
        "--clean_smiles",
        "-c",
        action="store_true",
        help="Pass SMILES through a set of standardizations. Default False.",
        dest="clean_smiles",
    )

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format="%(message)s")

    # elif args.subcommand == "run":
    logging.info("Beginning run...")
    start = time.time()
    run(
        ids_path=args.ids_path,
        id_type=args.ids_type,
        output_path=args.output_path,
        get_related_chebis=args.complete_rgroups,
        clean_smiles=args.clean_smiles,
    )
    end = time.time()
    logging.info("Run finished.\nTotal time: {}".format(end - start))
