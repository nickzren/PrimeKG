import os
import pandas as pd

# Constants
OUTPUT_DIR = 'data/output/drugbank/'

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load the drug bonds and drug vocabulary data
df_bonds = pd.read_csv(os.path.join(OUTPUT_DIR, 'drug_bonds.tsv'), sep='\t')
df_vocab = pd.read_csv(os.path.join(OUTPUT_DIR, 'drugbank_vocabulary.tsv'), sep='\t')
gene_vocab = pd.read_csv('data/input/vocab/gene_map.csv', delimiter='\t')

# Create dictionaries for mapping
up2ncbi = gene_vocab.set_index('UniProt ID(supplied by UniProt)')['NCBI Gene ID(supplied by NCBI)'].to_dict()
dbid2name = df_vocab.set_index('drugbank_id')['name'].to_dict()

def add_col(df, dct, source, sink, dtype=None):
    df[sink] = df[source].map(dct).fillna(float("nan"))
    if dtype:
        df[sink] = df[sink].astype(dtype)
    return df

df_bonds = add_col(df_bonds, up2ncbi, 'uniprot_id', 'NCBIGeneID', 'Int64')

def drugbank2edges(df, relation):
    db = []
    for _, data in df.iterrows():
        drugs = data['Drug IDs'].split('; ')
        for drug in drugs:
            info = {
                'DrugBank': drug,
                'relation': relation,
                'NCBIGeneID': data['NCBIGeneID'],
                'UniProtName': data['name'],
                'UniProtID': data['uniprot_id'],
            }
            db.append(info)
    return db

db = []
for bond_type in ['carrier', 'enzyme', 'target', 'transporter']:
    df_bond_type = df_bonds[df_bonds['bond_type'].str.lower().str.contains(bond_type)]
    db.extend(drugbank2edges(df_bond_type, bond_type))

df_prot_drug = pd.DataFrame(db)
df_prot_drug = add_col(df_prot_drug, dbid2name, 'DrugBank', 'DrugBankName')

df_prot_drug.to_csv(os.path.join(OUTPUT_DIR, 'drug_protein.tsv'), sep='\t', index=False)

print("Data saved to data/output/drugbank/drug_protein.tsv")