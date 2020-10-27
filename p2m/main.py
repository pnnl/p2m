#%% Main with argparse options and execution

from rdkit import RDLogger
import argparse
import sys
import os
import time
from . import annotation, data_download

def main():
    parser = argparse.ArgumentParser(description="Identify metabolites (substrates and products) of proteins "
                                 "from UniProt identifiers or EC numbers. "
                                 "Requires internet access for external database queries. "
                                 "\nTo set up local databases for your first run, be sure to run the 'download' command.")
    subparsers = parser.add_subparsers(title='Subcommands', dest='subcommand', required=True)

    subparser_download = subparsers.add_parser('download', help='Download databases. Required for first time use.')

    subparser_download.add_argument('--download_path', type=str,
                       help='Download required databases to the --path location.',
                       dest='download_path',
                       required=True)

    subparser_run = subparsers.add_parser('run', help='Run P2M on a set of protein identifiers.')
    subparser_run.add_argument('--ids_path', type=str,
                       help='Path to UniProt or EC identifiers. This should be a list of identifiers separated by newlines.',
                       dest='ids_path',
                       required=True)
    subparser_run.add_argument('--ids_type', type=str,
                       help='Type of identifiers. One of: UniProt, EC',
                       dest='ids_type',
                       required=True)
    subparser_run.add_argument('--database_path', type=str,
                        help='Path to local Rhea and ChEBI databases.',
                        dest='database_path',
                        required=True)
    subparser_run.add_argument('--output_path', type=str,
                        help='Path to desired output folder.',
                        dest='output_path',
                        required=True)
    subparser_run.add_argument('--external_uniprot', action='store_true',
                        help='set to "True" to perform external mapping of UniProt identifiers without local Rhea crossreferences to the UniProt database. Default False.',
                        dest='external_uniprot')
    subparser_run.add_argument('--external_rhea', action='store_true',
                        help='set to "True" to attempt to include preliminary Rhea reactions by performing external Rhea database queries. Default False.',
                        dest='external_rhea')
    subparser_run.add_argument('--complete_rgroups', action='store_true',
                        help='set to "True" to externally query the ChEBI database for substructure searches of compounds with R-groups. Default False.',
                        dest='complete_rgroups')
    subparser_run.add_argument('--clean_smiles', action='store_true',
                        help='set to "True" to pass SMILES through a set of standardizations. Default False.',
                        dest='clean_smiles')


    args = parser.parse_args()
    print(args)
    if 'download_path' in args:
        print("Beginning downloads...")
        data_download.download_rhea_data(args.download_path)
        print("Complete. Data directory is located at: {}".format(args.download_path))
        sys.exit()

    print("Beginning run...")
    start = time.time()
    results = run(ids_path = args.ids_path,
                  id_type = args.ids_type,
                  data_path = args.database_path,
                  output_path = args.output_path,
                  external_rhea = args.external_rhea,
                  get_related_chebis = args.complete_rgroups,
                  query_external_uniprot = args.external_uniprot,
                  clean_smiles = args.clean_smiles)
    end = time.time()
    print("Run finished. Total time: {}".format(end-start))
    return(results)

def run(ids_path,
        id_type,
        data_path,
        output_path,
        external_rhea=False,
        get_related_chebis=False,
        query_external_uniprot=False,
        clean_smiles=False):
    """
    Perform workflow
    """
    tsv_path = os.path.join(data_path, 'tsv')
    tsv_files = ['rhea2ec.tsv',
                  'rhea-directions.tsv',
                  'rhea2uniprot.tsv',
                  'chebiId_name.tsv']
    rxn_path = os.path.join(data_path, 'rxn')
    tsv_file_paths = [os.path.join(tsv_path, file) for file in tsv_files]
    paths = dict(zip(tsv_files, tsv_file_paths))
    paths['rxn'] = rxn_path
    paths['data'] = data_path
    paths['chebi'] = os.path.join(data_path, 'ChEBI_complete.sdf')
    paths['chebi_smiles'] = os.path.join(data_path, 'ChEBI_complete_smiles.txt')

    RDLogger.DisableLog('rdApp.*')
    
    with open(ids_path) as f:
        ids = f.read().splitlines()
    
    annot = annotation.Annotation(ids, id_type, paths)
    print("Mapping protein identifiers to Rhea identifiers...")
    
    annot.map_rhea(query_external_uniprot)
    if external_rhea == True:
        print("Externally querying Rhea database for identifiers not found locally...")
        annot.add_external_rhea()
    annot.add_internal_rhea()
    print("Obtaining metabolites from local Rhea reactions...")
    annot.get_rxn_products()
    if get_related_chebis == True:
        print("Externally querying ChEBI database to add unspecified R-group containing compounds...")
        annot.star_to_smiles()
    if clean_smiles == True:
        print("Standardizing SMILES strings...")
        annot.clean_smiles()
    print("Exporting output to: {}".format(output_path))
    if os.path.exists(output_path) == False:
        os.mkdir(output_path)
    annot.export_smiles(output_path)
    return(annot)





