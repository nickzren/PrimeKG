import os
import pandas as pd

# Constants
OUTPUT_DIR = 'data/output/drugcentral'
INPUT_DIR = 'data/input/drugcentral/'

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Import drug data (from drugbank; identifier: DB ID)
# drug_list = pd.read_csv('drug_nodes_in_kgv3.csv')
# print(drug_list.head())
# print(len(drug_list))  # 7928

# Drug bank All data
db_all = pd.read_csv('data/output/drugbank/drugs.tsv', sep='\t')
print(db_all.head())
db_all = db_all.rename(columns={'drugbank_id': 'id', 'cas_number': 'CAS_ID', 'type': 'drug_type'})

# MAPPING TABLE: DB ID--CAS ID
db_mapping = db_all[['id', 'CAS_ID']]
print(db_mapping.head())
print(len(db_mapping))  # 14315

# Drug type
db_drug_type = db_all[['id', 'CAS_ID', 'drug_type']]
print(db_drug_type.head())
print(len(db_drug_type))  # 14315

# Drug Central Features

# Description (MW, TPSA, CLOGP)
dc_description = pd.read_csv(f'{INPUT_DIR}/structures.tsv', sep='\t')
print(dc_description.head())
# MW
dc_mw = dc_description[['cd_molweight', 'cas_reg_no', 'name', 'id']]
# TPSA
dc_tpsa = dc_description[['tpsa', 'cas_reg_no', 'name', 'id']]
# CLOGP
dc_clogp = dc_description[['clogp', 'cas_reg_no', 'name', 'id']]

# MAPPING TABLE: Structure ID--CAS ID
dc_mapping = dc_description[['id', 'cas_reg_no', 'name']].rename(columns={'id': 'Structure_ID', 'cas_reg_no': 'CAS_ID'})
print(len(dc_mapping))  # 4531

# MW
dc_mw = dc_mw[dc_mw['cd_molweight'] != 'NULL'].rename(columns={'cd_molweight': 'MW', 'cas_reg_no': 'CAS_ID', 'id': 'Structure_ID'})[['MW', 'Structure_ID']]
print(len(dc_mw))

# TPSA
dc_tpsa = dc_tpsa[dc_tpsa['tpsa'] != 'NULL'].rename(columns={'tpsa': 'TPSA', 'cas_reg_no': 'CAS_ID', 'id': 'Structure_ID'})[['TPSA', 'Structure_ID']]

# CLOGP
dc_clogp = dc_clogp[dc_clogp['clogp'] != 'NULL'].rename(columns={'clogp': 'CLOGP', 'cas_reg_no': 'CAS_ID', 'id': 'Structure_ID'})[['CLOGP', 'Structure_ID']]

# Structure type
dc_str_type = pd.read_csv(f'{INPUT_DIR}/structure_type.tsv', sep='\t')
print(dc_str_type.head())
dc_str_type = dc_str_type[['struct_id', 'type']].rename(columns={'struct_id': 'Structure_ID', 'type': 'STRUCTURE_TYPE'})
print(dc_str_type.head())

# Map structure id to CAS ID
dc_feature = dc_mapping.merge(dc_str_type, on='Structure_ID', how='left')
dc_feature = dc_feature.merge(dc_mw, on='Structure_ID', how='left')
dc_feature = dc_feature.merge(dc_tpsa, on='Structure_ID', how='left')
dc_feature = dc_feature.merge(dc_clogp, on='Structure_ID', how='left')
print(dc_feature.head())
print(len(dc_feature))

# Map CAS ID to DB ID
db_dc = db_mapping.merge(dc_feature, on='CAS_ID', how='left')
print(db_dc.head())
print(len(db_dc))

# Join drug_list with db_dc
# result = drug_list.merge(db_dc, on='id', how='left')
# print(result.head())

# Pka
# dc_pka = pd.read_csv('data/input/drugcentral/pka.csv')
# print(dc_pka.head())
# print(len(dc_pka))  # 7301
# print(len(dc_pka['struct_id'].unique()))  # 3344

# Save dc_feature to TSV file
dc_feature.to_csv(f'{OUTPUT_DIR}/dc_features.tsv', sep='\t', index=False)
print(f"Data saved to {OUTPUT_DIR}/dc_features.tsv")
