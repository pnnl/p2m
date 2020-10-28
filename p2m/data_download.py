import os
import pandas as pd
from rdkit import Chem
import urllib.request
import tarfile
import gzip
import shutil
from . import utils as u

def download_rhea_data(data_path):
    if os.path.exists(data_path):
        raise Exception("{} already exists. Please specify a new directory to download data to.".format(data_path))
    else:
        os.mkdir(data_path)
        tsv_path = os.path.join(data_path, 'tsv')
        os.mkdir(tsv_path)
    tsv_url = 'ftp://ftp.ebi.ac.uk/pub/databases/rhea/tsv/'
    rxn_url = 'ftp://ftp.expasy.org/databases/rhea/ctfiles/rhea-rxn.tar.gz'
    chebi_url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete.sdf.gz'
    tsv_files = ['rhea2ec.tsv',
             'rhea-directions.tsv',
             'rhea2uniprot.tsv',
             'chebiId_name.tsv']
    print("Downloading Rhea .tsv files...")
    [urllib.request.urlretrieve(os.path.join(tsv_url, file), os.path.join(tsv_path, file))[0] for file in tsv_files]
    
    print("Downloading compressed Rhea reactions...")
    rhea_gz = os.path.join(data_path, os.path.basename(rxn_url))
    urllib.request.urlretrieve(rxn_url, rhea_gz)
    
    print("Extracting Rhea reaction files...")
    tar_file = tarfile.open(rhea_gz, "r:gz")
    tar_file.extractall(path=data_path)
    tar_file.close()
    
    
    print("Downloading ChEBI complete database...")
    chebi_gz = os.path.join(data_path, os.path.basename(chebi_url))
    urllib.request.urlretrieve(chebi_url, chebi_gz)
    
    print("Extracting ChEBI database...")
    sdf_path = os.path.join(data_path, 'ChEBI_complete.sdf')
    with gzip.open(chebi_gz, 'rb') as chebi_in:
        with open(sdf_path, 'wb') as chebi_out:
            shutil.copyfileobj(chebi_in, chebi_out)
    
    print("Generating ChEBI ID to SMILES mapping file...")
    chebi_out_path = os.path.join(data_path, 'ChEBI_complete_smiles.txt') 
    suppl = Chem.SDMolSupplier(sdf_path)
    all_smiles = [u.get_chebi_info(x) for x in suppl if type(x) == Chem.rdchem.Mol]
    df = pd.DataFrame(all_smiles)
    df.columns = ['ChEBI_ID', 'Name', 'SMILES']
    df.to_csv(chebi_out_path, sep='\t', index=False)
    
    print("Cleaning up...")
    os.remove(rhea_gz)
    os.remove(chebi_gz)
    os.remove(sdf_path)

