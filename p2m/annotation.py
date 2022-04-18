from . import utils as u
import os
import pandas as pd

class Annotation:
    """Class for obtaining metabolites associated with UniProt or EC identifiers cross-referencing to the Rhea database."""
    def __init__(self, ids, ids_type, paths):
        self.ids = {}
        self.ids = ids
        self.ids_type = ids_type
        self.paths = paths
        self.rhea_ids = []
        self.rxns = []
    
    def map_uniprot(self, external_uniprot):
        """Map UniProt IDs to local Rhea reaction identifiers with optional external database querying"""
        rhea2uniprot = pd.read_csv(self.paths['rhea2uniprot.tsv'], sep='\t', dtype=str)
        found_uni_rheas, self.missing_uni_ids = u.get_master_rhea(rhea2uniprot, self.ids)
        self.rhea_ids.extend(found_uni_rheas)
        if external_uniprot==True:
            self.external_uniprot_rhea = u.query_external_uniprot(self.missing_uni_ids)
            self.rhea_ids.extend(self.external_uniprot_rhea)
    
    def map_ec(self):
        """Map EC IDs to local Rhea reaction identifiers"""
        rhea2ec = pd.read_csv(self.paths['rhea2ec.tsv'], sep='\t', dtype=str)
        found_ec_rheas, self.missing_ec_ids = u.get_master_rhea(rhea2ec, self.ids)
        self.rhea_ids.extend(found_ec_rheas)
        
    def map_rhea(self, external_uniprot):
        """Map ID's (UniProt or EC) to Rhea reaction identifiers."""
        if self.ids_type == 'UniProt':
            self.map_uniprot(external_uniprot)
        elif self.ids_type=='EC':
            self.map_ec()
        self.unique_rhea_ids = list(set(self.rhea_ids))   
        self.rhea_directions = pd.read_csv(self.paths['rhea-directions.tsv'], sep='\t', dtype=str)
        self.rhea_rxn_ids = u.get_rxn_ids(self.unique_rhea_ids, self.rhea_directions)
                
    def add_external_rhea(self):
        """
        Externally query the Rhea database for rxn files from Rhea reaction ids.
        Note that this guesses external directional ID's and adds both directional reactions for downstream id matching
        """
        guess_lr, guess_rl = [[str(int(x) + j) for x in self.rhea_rxn_ids['external']] for j in [1, 2]]
        all_externals = self.rhea_rxn_ids['external'] + guess_lr + guess_rl
        dir_externals = self.rhea_directions.isin(all_externals)
        invalid_externals = list(self.rhea_directions[dir_externals == True].stack().values)
        searchable_rxns = list(set(all_externals) - set(invalid_externals))
        self.external_rhea_rxns = u.external_rhea(searchable_rxns)
        self.rxns.extend(self.external_rhea_rxns)
            
    def add_internal_rhea(self):
        """Load rxn files based on internal rhea ids"""
        rxn_paths = [os.path.join(self.paths['rxn'], rhea_id) + '.rxn' for rhea_id in self.rhea_rxn_ids['local']]
        self.internal_rhea_rxns = [u.convert_rxn(path) for path in rxn_paths]
        self.rxns.extend(self.internal_rhea_rxns)

    def get_rxn_products(self):
        """Get product SMILES strings from .rxn files"""
        all_smiles = [u.rxn_to_smiles(rxn) for rxn in self.rxns]
        smiles_valid = [smiles for smiles in all_smiles if type(smiles) == dict]
        flat_smiles = {}
        for entry in smiles_valid:
            for chebi, smiles in entry.items():
                flat_smiles[chebi] = smiles
        smiles_df = pd.DataFrame.from_dict(flat_smiles, orient='index', columns=['smiles'])
        smiles_df.reset_index(inplace=True)
        smiles_df.columns = ['ChEBI_ID', 'SMILES']
        smiles_named = u.get_chebi_names(self.paths['chebiId_name.tsv'], smiles_df)
        self.smiles = smiles_named.drop_duplicates(subset='ChEBI_ID', ignore_index=True)
        self.smiles_stars = smiles_named[smiles_named.SMILES.str.contains("\*")]   
        self.smiles_nostars = smiles_named[~smiles_named.SMILES.str.contains("\*")]
    
    def star_to_smiles(self):
        """Externally query ChEBI database for related compounds"""
        related_chebis = [u.chebi_rels(chebi) for chebi in list(self.smiles_stars['ChEBI_ID']) if "CHEBI" in chebi]
        flattened_related = ["CHEBI:" + str(i) for s in related_chebis for i in s]
        matched_chebis = u.get_chebi_smiles(pd.read_csv(self.paths['chebi_smiles'], sep='\t'), flattened_related)
        all_smiles = matched_chebis.drop_duplicates(subset='ChEBI_ID', ignore_index=True)
        nostars = pd.concat([self.smiles_nostars, all_smiles[~all_smiles.SMILES.str.contains("\*")]])
        self.all_nostars = nostars.drop_duplicates(subset='ChEBI_ID', ignore_index=True)   
    
    def clean_smiles(self):
        """
        Pass all SMILES strings in the dataframe through standardization steps 
        that desalt, uncharge, and generate a canonical tautomer for each
        SMILES entry in the dataframe.
        """
        self.clean_smiles_nostars = u.clean_smiles_df(self.smiles_nostars)
        if hasattr(self, 'all_nostars'):
            self.clean_all_nostars = u.clean_smiles_df(self.all_nostars)
    
    def export_smiles(self, export_path):
        """Export SMILES dataframes with ChEBI identifiers to appropriate locationsin the output folder."""
        self.smiles_stars.to_csv(os.path.join(export_path, 'smiles_stars.tsv'), index=False, sep='\t')
        if hasattr(self, 'clean_smiles_nostars'):
            self.clean_smiles_nostars.to_csv(os.path.join(export_path, 'clean_smiles_nostars.tsv'), index=False, sep='\t')
            if hasattr(self, 'clean_all_nostars'):
                self.clean_all_nostars.to_csv(os.path.join(export_path, 'clean_all_nostars.tsv'), index=False, sep='\t')
        else:
            self.smiles_nostars.to_csv(os.path.join(export_path, 'smiles_nostars.tsv'), index=False, sep='\t')
            if hasattr(self, 'all_nostars'):
                self.all_nostars.to_csv(os.path.join(export_path, 'smiles_all_nostars.tsv'), index=False, sep='\t')
     
                                              
