import numpy as np
import pandas as pd

# Constants for file paths
MONDO_TERMS_FILE = 'data/output/mondo/mondo_terms.csv'
MONDO_XREF_FILE = 'data/output/mondo/mondo_references.csv'
UMLS_FILE = 'data/output/umls/umls.tsv'
OUTPUT_FILE = 'data/output/vocab/umls_mondo.csv'

# Valid source-ontology pairs for indirect mapping
VALID_PAIRS = [
    ('OMIM', 'OMIM'),
    # ('OMIM', 'OMIMPS'),
    # ('NCI_NICHD', 'NCIT'),
    # ('NCI_FDA', 'NCIT'),
    # ('NCI_CTRP', 'NCIT'),
    ('NCI', 'NCIT'),
    ('MSH', 'MESH'),
    ('MDR', 'MedDRA'),
    # ('ICD10CM', 'ICD10'),
    ('ICD10', 'ICD10'),
    ('SNOMEDCT_US', 'SCTID')
]

def load_data(file_path, columns=None, sep=',', dtype=None):
    """Load data from a CSV file."""
    return pd.read_csv(file_path, usecols=columns, sep=sep, dtype=dtype)

def get_direct_mapping(df_mondo_xref):
    """Get direct mapping from MONDO references to UMLS."""
    return df_mondo_xref.query('ontology == "UMLS"')[['ontology_id', 'mondo_id']].rename(columns={'ontology_id': 'umls_id'})

def get_umls_joined(df_umls):
    """Get UMLS data with selected columns and drop duplicates."""
    return df_umls[['cui', 'source_cui', 'source_descriptor_dui', 'source', 'source_code']].drop_duplicates()

def get_indirect_mapping(df_umls_join, df_mondo_xref, valid_pairs):
    """Get indirect mapping from UMLS to MONDO using valid source-ontology pairs."""
    map1 = pd.merge(df_umls_join, df_mondo_xref, how='inner', left_on='source_cui', right_on='ontology_id')
    map2 = pd.merge(df_umls_join, df_mondo_xref, how='inner', left_on='source_descriptor_dui', right_on='ontology_id')
    map3 = pd.merge(df_umls_join, df_mondo_xref, how='inner', left_on='source_code', right_on='ontology_id')
    
    map_all = pd.concat([map1, map2, map3]).drop_duplicates()

    valid_mappings = []
    for source, ontology in valid_pairs:
        valid_mappings.append(map_all.query('source == @source and ontology == @ontology'))
    
    return pd.concat(valid_mappings).rename(columns={'cui': 'umls_id'})

def combine_mappings(map_direct, map_indirect):
    """Combine direct and indirect mappings and drop duplicates."""
    return pd.concat([map_direct, map_indirect[['umls_id', 'mondo_id']]]).drop_duplicates()

def main():
    # Load data
    df_mondo_terms = load_data(MONDO_TERMS_FILE, columns=['name', 'id'])
    df_mondo_xref = load_data(MONDO_XREF_FILE)
    dtype_spec = {'cui': str, 'source_cui': str, 'source_descriptor_dui': str, 'source': str, 'source_code': str}
    df_umls = load_data(UMLS_FILE, sep='\t', dtype=dtype_spec)

    # Get mappings
    map_direct = get_direct_mapping(df_mondo_xref)
    df_umls_join = get_umls_joined(df_umls)
    map_indirect = get_indirect_mapping(df_umls_join, df_mondo_xref, VALID_PAIRS)

    # Combine mappings
    df_umls_mondo = combine_mappings(map_direct, map_indirect)

    # Save the result to a CSV file
    df_umls_mondo.to_csv(OUTPUT_FILE, index=False)
    print(f"Processing complete. Saved to {OUTPUT_FILE}")

if __name__ == '__main__':
    main()
