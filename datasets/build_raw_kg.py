import os
import pandas as pd
import numpy as np
import igraph as ig

# Constants
INPUT = 'data/input/'
OUTPUT = 'data/output/'
SAVE_PATH = OUTPUT + 'kg/'

def assert_dtypes(df): 
    all_string = True
    for i, x in enumerate(df.dtypes.values): 
        if x != np.dtype('O'): 
            all_string = False
            print(df.columns[i], x)
    if not all_string: assert False

def clean_edges(df): 
    df = df.get(['relation', 'display_relation', 'x_id','x_type', 'x_name', 'x_source','y_id','y_type', 'y_name', 'y_source'])
    df = df.dropna()
    df = df.drop_duplicates()
    df = df.query('not ((x_id == y_id) and (x_type == y_type) and (x_source == y_source) and (x_name == y_name))')
    return df

# Load DataFrames
df_drugbank = pd.read_csv(OUTPUT + 'drugbank/drug_protein.tsv', sep='\t', low_memory=False)
df_drugbank = df_drugbank.get(['DrugBank', 'relation', 'NCBIGeneID','DrugBankName']).dropna()
df_drugbank = df_drugbank.astype({'NCBIGeneID':int}).astype({'NCBIGeneID':str})
assert_dtypes(df_drugbank)

df_disgenet = pd.read_csv(INPUT + 'disgenet/curated_gene_disease_associations.tsv', sep='\t', low_memory=False)
df_disgenet = df_disgenet.astype({'geneId':int}).astype({'geneId':str})

df_mondo_terms = pd.read_csv(OUTPUT + 'mondo/mondo_terms.csv', low_memory=False)
df_mondo_terms = df_mondo_terms[pd.to_numeric(df_mondo_terms['id'], errors='coerce').notnull()]
df_mondo_terms = df_mondo_terms.astype({'id': int}).astype({'id': str})

df_mondo_xref = pd.read_csv(OUTPUT + 'mondo/mondo_references.csv', low_memory=False)
df_mondo_xref = df_mondo_xref.astype({'mondo_id':int}).astype({'mondo_id':str})
assert_dtypes(df_mondo_xref)

df_mondo_parents = pd.read_csv(OUTPUT + 'mondo/mondo_parents.csv', low_memory=False)
df_mondo_parents = df_mondo_parents[pd.to_numeric(df_mondo_parents['parent'], errors='coerce').notnull()]
df_mondo_parents = df_mondo_parents[pd.to_numeric(df_mondo_parents['child'], errors='coerce').notnull()]
df_mondo_parents = df_mondo_parents.astype({'parent':int}).astype({'parent':str})
df_mondo_parents = df_mondo_parents.astype({'child':int}).astype({'child':str})
assert_dtypes(df_mondo_parents)

df_drug_central = pd.read_csv(OUTPUT + 'drugcentral/drug_disease.csv', low_memory=False)
df_drug_central = df_drug_central.get(['cas_reg_no','relationship_name', 'umls_cui']) # 'concept_id', 'concept_name', 'snomed_conceptid'
df_drug_central = df_drug_central.query('not @df_drug_central.cas_reg_no.isna()')
df_drug_central = df_drug_central.query('not @df_drug_central.umls_cui.isna()')
assert_dtypes(df_drug_central)

df_ddi = pd.read_csv(OUTPUT + 'drugbank/drug_drug.tsv', sep='\t', low_memory=False)
assert_dtypes(df_ddi)

df_hp_terms = pd.read_csv(OUTPUT + 'hpo/hp_terms.csv', low_memory=False)
df_hp_terms = df_hp_terms.astype({'id':int}).astype({'id':str})

df_hp_xref = pd.read_csv(OUTPUT + 'hpo/hp_references.csv', low_memory=False)
df_hp_xref = df_hp_xref.astype({'hp_id':int}).astype({'hp_id':str})

df_hp_parents = pd.read_csv(OUTPUT + 'hpo/hp_parents.csv', low_memory=False)
df_hp_parents = df_hp_parents.astype({'parent':int}).astype({'parent':str})
df_hp_parents = df_hp_parents.astype({'child':int}).astype({'child':str})
assert_dtypes(df_hp_parents)

df_hpoa_pos = pd.read_csv(OUTPUT + 'hpo/disease_phenotype_pos.csv', low_memory=False)
df_hpoa_pos = df_hpoa_pos.astype({'hp_id':int}).astype({'hp_id':str})
df_hpoa_pos = df_hpoa_pos.astype({'disease_ontology_id':int}).astype({'disease_ontology_id':str})
assert_dtypes(df_hpoa_pos)

df_hpoa_neg = pd.read_csv(OUTPUT + 'hpo/disease_phenotype_neg.csv', low_memory=False)
df_hpoa_neg = df_hpoa_neg.astype({'hp_id':int}).astype({'hp_id':str})
df_hpoa_neg = df_hpoa_neg.astype({'disease_ontology_id':int}).astype({'disease_ontology_id':str})
assert_dtypes(df_hpoa_neg)

df_sider = pd.read_csv(OUTPUT + 'sider/sider.csv', low_memory=False)
assert_dtypes(df_sider)

df_go_terms = pd.read_csv(OUTPUT + 'go/go_terms_info.csv', low_memory=False)
df_go_terms = df_go_terms.astype({'go_term_id':int}).astype({'go_term_id':str})
assert_dtypes(df_go_terms)

df_go_edges = pd.read_csv(OUTPUT + 'go/go_terms_relations.csv', low_memory=False)
df_go_edges = df_go_edges.astype({'x':int}).astype({'x':str})
df_go_edges = df_go_edges.astype({'y':int}).astype({'y':str})
assert_dtypes(df_go_edges)

df_gene2go = pd.read_csv(OUTPUT + 'ncbigene/protein_go_associations.csv', low_memory=False)
df_gene2go = df_gene2go.astype({'ncbi_gene_id':int}).astype({'ncbi_gene_id':str})
df_gene2go = df_gene2go.astype({'go_term_id':int}).astype({'go_term_id':str})
assert_dtypes(df_gene2go)

df_exposures = pd.read_csv(OUTPUT + 'ctd/exposure_data.csv', low_memory=False)
df_exposures = df_exposures.get(['exposurestressorname', 'exposurestressorid',
                  'exposuremarker', 'exposuremarkerid',
                  'diseasename', 'diseaseid',
                  'phenotypename', 'phenotypeid'])
assert_dtypes(df_exposures)

df_uberon_terms = pd.read_csv(OUTPUT + 'uberon/uberon_terms.csv', low_memory=False)
df_uberon_terms = df_uberon_terms.astype({'id':int}).astype({'id':str})
assert_dtypes(df_uberon_terms)

df_uberon_is_a = pd.read_csv(OUTPUT + 'uberon/uberon_is_a.csv', low_memory=False)
df_uberon_is_a = df_uberon_is_a.astype({'id':int}).astype({'id':str})
df_uberon_is_a = df_uberon_is_a.astype({'is_a':int}).astype({'is_a':str})
assert_dtypes(df_uberon_is_a)

df_uberon_rels = pd.read_csv(OUTPUT + 'uberon/uberon_rels.csv', low_memory=False)
df_uberon_rels = df_uberon_rels.astype({'id':int}).astype({'id':str})
df_uberon_rels = df_uberon_rels.astype({'relation_id':int}).astype({'relation_id':str})
assert_dtypes(df_uberon_rels)

df_bgee = pd.read_csv(OUTPUT + 'bgee/anatomy_gene.csv', low_memory=False)
df_bgee = df_bgee.astype({'expression_rank':int}).astype({'expression_rank':str})
df_bgee = df_bgee.astype({'anatomy_id':int}).astype({'anatomy_id':str})
assert_dtypes(df_bgee)

df_reactome_terms = pd.read_csv(OUTPUT + 'reactome/reactome_terms.csv', low_memory=False)
assert_dtypes(df_reactome_terms)

df_reactome_rels = pd.read_csv(OUTPUT + 'reactome/reactome_relations.csv', low_memory=False)
assert_dtypes(df_reactome_rels)

df_reactome_ncbi = pd.read_csv(OUTPUT + 'reactome/reactome_ncbi.csv', low_memory=False)
df_reactome_ncbi = df_reactome_ncbi[df_reactome_ncbi.ncbi_id.str.isnumeric()]
assert_dtypes(df_reactome_ncbi)

df_umls_mondo = pd.read_csv(OUTPUT + 'vocab/umls_mondo.csv', low_memory=False)
df_umls_mondo = df_umls_mondo.astype({'mondo_id':int}).astype({'mondo_id':str})
assert_dtypes(df_umls_mondo)

df_prot_names = pd.read_csv(INPUT + 'vocab/gene_names.csv', low_memory=False, sep='\t')
df_prot_names = df_prot_names.rename(columns={'NCBI Gene ID(supplied by NCBI)':'ncbi_id', 'NCBI Gene ID':'ncbi_id2', 'Approved symbol':'symbol', 'Approved name':'name'})
df_prot_names = df_prot_names.get(['ncbi_id', 'symbol']).dropna()
df_prot_names = df_prot_names.astype({'ncbi_id':int}).astype({'ncbi_id':str})
assert_dtypes(df_prot_names)

db_vocab = pd.read_csv(OUTPUT+'drugbank/drugbank_vocabulary.tsv', sep='\t', low_memory=False)
assert_dtypes(db_vocab)

df_db_atc = pd.read_csv(OUTPUT+'drugbank/drugbank_atc_codes.tsv', sep='\t', low_memory=False).get(['atc_code','parent_key'])
assert_dtypes(df_db_atc)

### Drug protein interactions (DrugBank)
df_prot_drug = pd.merge(df_drugbank, df_prot_names, 'left', left_on='NCBIGeneID', right_on='ncbi_id')

df_prot_drug = df_prot_drug.rename(columns={'DrugBank':'x_id', 'NCBIGeneID':'y_id', 'DrugBankName':'x_name', 'symbol':'y_name'})
df_prot_drug['x_type'] = 'drug'
df_prot_drug['x_source'] = 'DrugBank'
df_prot_drug['y_type'] = 'gene/protein'
df_prot_drug['y_source'] = 'NCBI'
df_prot_drug['display_relation'] = df_prot_drug.get('relation').values
df_prot_drug['relation'] = 'drug_protein' # combine targets, carrier, enzyme and transporter
df_prot_drug = clean_edges(df_prot_drug)
df_prot_drug.head(1)

# Drug disease interactions (DiseaseCentral)
df_drug_dis = pd.merge(df_drug_central, db_vocab, 'left', left_on='cas_reg_no', right_on='cas_number')
df_drug_dis = pd.merge(df_drug_dis, df_umls_mondo, 'inner', left_on='umls_cui', right_on='umls_id')
df_drug_dis = pd.merge(df_drug_dis, df_mondo_terms, 'left', left_on='mondo_id', right_on='id')

df_drug_dis = df_drug_dis.get(['relationship_name', 'drugbank_id', 'name_x', 'mondo_id', 'name_y'])
df_drug_dis = df_drug_dis.dropna().drop_duplicates()

df_drug_dis = df_drug_dis.rename(columns={'drugbank_id':'x_id', 'mondo_id':'y_id', 'name_x':'x_name', 'name_y':'y_name', 'relationship_name':'relation'})
df_drug_dis['x_type'] = 'drug'
df_drug_dis['x_source'] = 'DrugBank'
df_drug_dis['y_type'] = 'disease'
df_drug_dis['y_source'] = 'MONDO'
df_drug_dis['display_relation'] = df_drug_dis.get('relation').values
df_drug_dis = clean_edges(df_drug_dis)
print(df_drug_dis.head(1))

### Disease protein interactions (DisGenNet)

df_prot_dis1 = df_disgenet.query('diseaseType=="disease"')

df_prot_dis1 = pd.merge(df_prot_dis1, df_umls_mondo, 'inner', left_on='diseaseId', right_on='umls_id')
df_prot_dis1 = pd.merge(df_prot_dis1, df_mondo_terms, 'left', left_on='mondo_id', right_on='id')

df_prot_dis1 = df_prot_dis1.rename(columns={'geneId':'x_id', 'geneSymbol':'x_name', 'mondo_id':'y_id', 'name':'y_name'})
df_prot_dis1['x_type'] = 'gene/protein'
df_prot_dis1['x_source'] = 'NCBI'
df_prot_dis1['y_type'] = 'disease'
df_prot_dis1['y_source'] = 'MONDO'
df_prot_dis1['relation'] = 'disease_protein'
df_prot_dis1['display_relation'] = 'associated with'
df_prot_dis1 = clean_edges(df_prot_dis1)
print(df_prot_dis1.head(1))

# Disease disease interations

df_dis_dis1 = pd.merge(df_mondo_parents, df_mondo_terms, 'left', left_on='parent', right_on='id')
df_dis_dis1 = df_dis_dis1.rename(columns={'parent':'x_id', 'name':'x_name'})
df_dis_dis1 = pd.merge(df_dis_dis1, df_mondo_terms, 'left', left_on='child', right_on='id')
df_dis_dis1 = df_dis_dis1.rename(columns={'child':'y_id', 'name':'y_name'})
df_dis_dis1['x_type'] = 'disease'
df_dis_dis1['x_source'] = 'MONDO'
df_dis_dis1['y_type'] = 'disease'
df_dis_dis1['y_source'] = 'MONDO'
df_dis_dis1['relation'] = 'disease_disease'
df_dis_dis1['display_relation'] = 'parent-child'
df_dis_dis1 = clean_edges(df_dis_dis1)
print(df_dis_dis1.head(1))

# Drug drug interactions (DrugBank)

df_drug_drug = pd.merge(df_ddi, db_vocab, 'inner', left_on='drug1', right_on='drugbank_id')
df_drug_drug = df_drug_drug.rename(columns={'drug1':'x_id', 'name':'x_name'})
df_drug_drug = pd.merge(df_drug_drug.astype({'drug2':'str'}), db_vocab, 'inner', left_on='drug2', right_on='drugbank_id')
df_drug_drug = df_drug_drug.rename(columns={'drug2':'y_id', 'name':'y_name'})
df_drug_drug['x_type'] = 'drug'
df_drug_drug['x_source'] = 'DrugBank'
df_drug_drug['y_type'] = 'drug'
df_drug_drug['y_source'] = 'DrugBank'
df_drug_drug['relation'] = 'drug_drug'
df_drug_drug['display_relation'] = 'synergistic interaction'
df_drug_drug = clean_edges(df_drug_drug)
print(df_drug_drug.head(1))

# EFFECT / Phenotype 
### Effect protein interactions (DisGenNet)

df_prot_phe = df_disgenet.query('diseaseType=="phenotype"')

df_prot_phe = pd.merge(df_prot_phe, df_hp_xref, 'inner', left_on='diseaseId', right_on='ontology_id')
df_prot_phe = pd.merge(df_prot_phe, df_hp_terms, 'left', left_on='hp_id', right_on='id')

df_prot_phe = df_prot_phe.rename(columns={'geneId':'x_id', 'geneSymbol':'x_name', 'hp_id':'y_id', 'name':'y_name'})
df_prot_phe['x_type'] = 'gene/protein'
df_prot_phe['x_source'] = 'NCBI'
df_prot_phe['y_type'] = 'effect/phenotype'
df_prot_phe['y_source'] = 'HPO'
df_prot_phe['relation'] = 'phenotype_protein'
df_prot_phe['display_relation'] = 'associated with'
df_prot_phe = clean_edges(df_prot_phe)
print(df_prot_phe.head(1))

### Effect effect interactions (HPO)

df_phe_phe = pd.merge(df_hp_parents, df_hp_terms, 'left', left_on='parent', right_on='id')
df_phe_phe = df_phe_phe.rename(columns={'name':'parent_name'})
df_phe_phe = pd.merge(df_phe_phe, df_hp_terms, 'left', left_on='child', right_on='id')
df_phe_phe = df_phe_phe.rename(columns={'name':'child_name'})
df_phe_phe = df_phe_phe.get(['parent', 'child', 'parent_name', 'child_name'])

df_phe_phe = df_phe_phe.rename(columns={'parent':'x_id', 'child':'y_id', 'parent_name':'x_name', 'child_name':'y_name'})
df_phe_phe['x_type'] = 'effect/phenotype'
df_phe_phe['x_source'] = 'HPO'
df_phe_phe['y_type'] = 'effect/phenotype'
df_phe_phe['y_source'] = 'HPO'
df_phe_phe['relation'] = 'phenotype_phenotype'
df_phe_phe['display_relation'] = 'parent-child'
df_phe_phe = clean_edges(df_phe_phe)
print(df_phe_phe.head(1))

### Disease effect interactions (HPO-A)

df_dis_phe_pos1 = pd.merge(df_hpoa_pos, df_mondo_xref, 'left', left_on='disease_ontology_id', right_on='ontology_id')
df_dis_phe_pos1 = df_dis_phe_pos1.query('(disease_ontology==ontology) or (disease_ontology=="ORPHA" and ontology=="Orphanet")')
df_dis_phe_pos1 = pd.merge(df_dis_phe_pos1, df_hp_terms, 'left', left_on='hp_id', right_on='id').rename(columns={'name':'hp_name'})
df_dis_phe_pos1 = pd.merge(df_dis_phe_pos1, df_mondo_terms, 'left', left_on='mondo_id', right_on='id').rename(columns={'name':'mondo_name'})
df_dis_phe_pos1 = df_dis_phe_pos1.get(['mondo_id', 'mondo_name', 'hp_id', 'hp_name'])
df_dis_phe_pos1 = df_dis_phe_pos1.rename(columns={'mondo_id':'x_id', 'mondo_name':'x_name', 'hp_id': 'y_id', 'hp_name':'y_name'})
df_dis_phe_pos1.loc[:, 'x_source'] = 'MONDO'
df_dis_phe_pos1.loc[:, 'x_type'] = 'disease'
df_dis_phe_pos1.loc[:, 'y_source'] = 'HPO'
df_dis_phe_pos1.loc[:, 'y_type'] = 'effect/phenotype'
df_dis_phe_pos1.loc[:, 'relation'] = 'disease_phenotype_positive'
df_dis_phe_pos1.loc[:, 'display_relation'] = 'phenotype present'
df_dis_phe_pos1 = clean_edges(df_dis_phe_pos1)
print(df_dis_phe_pos1.head(1))

df_dis_phe_neg = pd.merge(df_hpoa_neg, df_mondo_xref, 'left', left_on='disease_ontology_id', right_on='ontology_id')
df_dis_phe_neg = df_dis_phe_neg.query('(disease_ontology==ontology) or (disease_ontology=="ORPHA" and ontology=="Orphanet")')
df_dis_phe_neg = pd.merge(df_dis_phe_neg, df_hp_terms, 'left', left_on='hp_id', right_on='id').rename(columns={'name':'hp_name'})
df_dis_phe_neg = pd.merge(df_dis_phe_neg, df_mondo_terms, 'left', left_on='mondo_id', right_on='id').rename(columns={'name':'mondo_name'})
df_dis_phe_neg = df_dis_phe_neg.get(['mondo_id', 'mondo_name', 'hp_id', 'hp_name'])
df_dis_phe_neg = df_dis_phe_neg.rename(columns={'mondo_id':'x_id', 'mondo_name':'x_name', 'hp_id': 'y_id', 'hp_name':'y_name'})
df_dis_phe_neg.loc[:, 'x_source'] = 'MONDO'
df_dis_phe_neg.loc[:, 'x_type'] = 'disease'
df_dis_phe_neg.loc[:, 'y_source'] = 'HPO'
df_dis_phe_neg.loc[:, 'y_type'] = 'effect/phenotype'
df_dis_phe_neg.loc[:, 'relation'] = 'disease_phenotype_negative'
df_dis_phe_neg.loc[:, 'display_relation'] = 'phenotype absent'
df_dis_phe_neg = clean_edges(df_dis_phe_neg)
print(df_dis_phe_neg.head(1))

### Remove phenotype nodes if they exist in MONDO

# phenotypes that are actually diseases in MONDO
# avoid duplicate nodes and convert them to disease relations
mondo_xref_hp_subset = df_mondo_xref.query('ontology=="HP"')
mondo_xref_hp_subset.loc[:, 'ontology_id'] = mondo_xref_hp_subset.get('ontology_id').astype(int).astype(str).values
hp_ids_r_mondo = pd.merge(mondo_xref_hp_subset, df_hp_terms, 'inner', left_on='ontology_id', right_on='id').get('ontology_id').values

def replace_hp_data_w_mondo(df, hp_id_col, drop_cols=[]): 
    cols = list(df.columns.values)
    cols.extend(['mondo_id', 'mondo_name'])
    [cols.remove(x) for x in drop_cols]
    df = pd.merge(df, mondo_xref_hp_subset, 'left', left_on=hp_id_col, right_on='ontology_id')
    df = pd.merge(df, df_mondo_terms, 'left', left_on='mondo_id', right_on='id')
    df = df.rename(columns={'name':'mondo_name'}).get(cols)
    return df

# HANDLE EFFECT EFFECT 
# PHE-PHE should be PHE-DIS if ONE PHE is in MONDO

df_dis_phe_x = df_phe_phe.query('x_id in @hp_ids_r_mondo and y_id not in @hp_ids_r_mondo')
df_dis_phe_x = replace_hp_data_w_mondo(df=df_dis_phe_x, hp_id_col='x_id', 
                                       drop_cols=[c for c in df_dis_phe_x.columns.values if 'x_' in c])
df_dis_phe_x = df_dis_phe_x.rename(columns={'mondo_id':'x_id', 'mondo_name':'x_name'})
df_dis_phe_x.loc[:, 'x_source'] = 'MONDO'
df_dis_phe_x.loc[:, 'x_type'] = 'disease'

df_dis_phe_y = df_phe_phe.query('y_id in @hp_ids_r_mondo and x_id not in @hp_ids_r_mondo')
df_dis_phe_y = replace_hp_data_w_mondo(df=df_dis_phe_y, hp_id_col='y_id',
                                       drop_cols=[c for c in df_dis_phe_y.columns.values if 'y_' in c])
df_dis_phe_y = df_dis_phe_y.rename(columns={'mondo_id':'y_id', 'mondo_name':'y_name'})
df_dis_phe_y.loc[:, 'y_source'] = 'MONDO'
df_dis_phe_y.loc[:, 'y_type'] = 'disease'

df_dis_phe_pos2 = pd.concat([df_dis_phe_x, df_dis_phe_y], ignore_index=True)
df_dis_phe_pos2['relation'] = 'disease_phenotype_positive'
df_dis_phe_pos2.loc[:, 'display_relation'] = 'phenotype present'
df_dis_phe_pos2 = clean_edges(df_dis_phe_pos2)

# PHE-PHE should be DIS-DIS if BOTH PHE are in MONDO

df_dis_dis2 = df_phe_phe.query('x_id in @hp_ids_r_mondo and y_id in @hp_ids_r_mondo')
df_dis_dis2 = replace_hp_data_w_mondo(df=df_dis_dis2, 
                                       hp_id_col='x_id', 
                                       drop_cols=[c for c in df_dis_dis2.columns.values if 'x_' in c])
df_dis_dis2 = df_dis_dis2.rename(columns={'mondo_id':'x_id', 'mondo_name':'x_name'})
df_dis_dis2 = replace_hp_data_w_mondo(df=df_dis_dis2, 
                                       hp_id_col='y_id', 
                                       drop_cols=[c for c in df_dis_dis2.columns.values if 'y_' in c])
df_dis_dis2 = df_dis_dis2.rename(columns={'mondo_id':'y_id', 'mondo_name':'y_name'})
df_dis_dis2.loc[:, 'x_source'] = 'MONDO'
df_dis_dis2.loc[:, 'x_type'] = 'disease'
df_dis_dis2.loc[:, 'y_source'] = 'MONDO'
df_dis_dis2.loc[:, 'y_type'] = 'disease'
df_dis_dis2.loc[:,'relation'] = 'disease_disease'
df_dis_dis2.loc[:,'display_relation'] = 'parent-child'
df_dis_dis2 = clean_edges(df_dis_dis2)

# drop relations in PHE PHE if either PHE is in MONDO
# phenotype phenotype should have no disease nodes
df_phe_phe = df_phe_phe.query('x_id not in @hp_ids_r_mondo and y_id not in @hp_ids_r_mondo')

# HANDLE PROTEIN EFFECT 

# if phenotype in MONDO make it protein-disease relations 
df_prot_dis2= df_prot_phe.query('y_id in @hp_ids_r_mondo')
df_prot_dis2 = replace_hp_data_w_mondo(df=df_prot_dis2, hp_id_col='y_id',
                                       drop_cols=[c for c in df_prot_dis2.columns.values if 'y_' in c])
df_prot_dis2 = df_prot_dis2.rename(columns={'mondo_id':'y_id', 'mondo_name':'y_name'})
df_prot_dis2.loc[:, 'y_source'] = 'MONDO'
df_prot_dis2.loc[:, 'y_type'] = 'disease'
df_prot_dis2.loc[:, 'relation'] = 'disease_protein'
df_prot_dis2.loc[:, 'display_relation'] = 'associated with'
df_prot_dis2 = clean_edges(df_prot_dis2)

# remove from protein-phenotype if phenotype in MONDO 
df_prot_phe = df_prot_phe.query('y_id not in @hp_ids_r_mondo')

# HANDLE DISEASE EFFECT 

# remove from protein-phenotype if phenotype in MONDO 
df_dis_phe_pos1 = df_dis_phe_pos1.query('y_id not in @hp_ids_r_mondo')

# NEGATIVE disease_phenotype should just be dropped because negative disease_disease doesn't make sense 
df_dis_phe_neg = df_dis_phe_neg.query('y_id not in @hp_ids_r_mondo')

# COMBINE DATAFRAMES 

df_prot_dis = pd.concat([df_prot_dis1, df_prot_dis2], ignore_index=True).drop_duplicates()
df_dis_dis = pd.concat([df_dis_dis1, df_dis_dis2], ignore_index=True).drop_duplicates()
df_dis_phe_pos = pd.concat([df_dis_phe_pos1, df_dis_phe_pos2], ignore_index=True).drop_duplicates()

### Drug effect interactions (SIDER)

df_drug_effect = pd.merge(df_sider, df_db_atc, 'left', left_on='atc', right_on='atc_code')
df_drug_effect = df_drug_effect.rename(columns={'parent_key':'DrugBank', 'UMLS_from_meddra':'UMLS'})
df_drug_effect = pd.merge(df_drug_effect, db_vocab, 'left', left_on='DrugBank', right_on='drugbank_id')
df_drug_effect = pd.merge(df_drug_effect, df_hp_xref, 'left', left_on='UMLS' , right_on='ontology_id')
df_drug_effect = pd.merge(df_drug_effect, df_hp_terms, 'left', left_on='hp_id' , right_on='id')
df_drug_effect = df_drug_effect.get(['drugbank_id','name_x','hp_id', 'name_y'])
df_drug_effect = df_drug_effect.dropna().drop_duplicates()

df_drug_effect = df_drug_effect.rename(columns={'drugbank_id':'x_id', 'name_x':'x_name', 'hp_id':'y_id', 'name_y':'y_name'})
df_drug_effect['x_type'] = 'drug'
df_drug_effect['x_source'] = 'DrugBank'
df_drug_effect['y_type'] = 'effect/phenotype'
df_drug_effect['y_source'] = 'HPO'
df_drug_effect['relation'] = 'drug_effect'
df_drug_effect['display_relation'] = 'side effect'
df_drug_effect = df_drug_effect.query('y_id not in @hp_ids_r_mondo')
df_drug_effect = clean_edges(df_drug_effect)
print(df_drug_effect.head(1))

## GO Terms
### Go terms interactions (GO)

bp = df_go_terms.query('go_term_type=="biological_process"')
df_bp_bp = pd.merge(df_go_edges, bp, 'inner', left_on='x', right_on='go_term_id')
df_bp_bp = df_bp_bp.rename(columns={'go_term_id':'x_id','go_term_name':'x_name','go_term_type':'x_type'})
df_bp_bp = pd.merge(df_bp_bp, bp, 'inner', left_on='y', right_on='go_term_id')
df_bp_bp = df_bp_bp.rename(columns={'go_term_id':'y_id','go_term_name':'y_name','go_term_type':'y_type'})
df_bp_bp['relation'] = 'bioprocess_bioprocess'
df_bp_bp['x_source'] = 'GO'
df_bp_bp['y_source'] = 'GO'
df_bp_bp['display_relation'] = 'parent-child'
df_bp_bp = clean_edges(df_bp_bp)
print(df_bp_bp.head(1))

mf = df_go_terms.query('go_term_type=="molecular_function"')
df_mf_mf = pd.merge(df_go_edges, mf, 'inner', left_on='x', right_on='go_term_id')
df_mf_mf = df_mf_mf.rename(columns={'go_term_id':'x_id','go_term_name':'x_name','go_term_type':'x_type'})
df_mf_mf = pd.merge(df_mf_mf, mf, 'inner', left_on='y', right_on='go_term_id')
df_mf_mf = df_mf_mf.rename(columns={'go_term_id':'y_id','go_term_name':'y_name','go_term_type':'y_type'})
df_mf_mf['relation'] = 'molfunc_molfunc'
df_mf_mf['display_relation'] = 'parent-child'
df_mf_mf['x_source'] = 'GO'
df_mf_mf['y_source'] = 'GO'
df_mf_mf = clean_edges(df_mf_mf)
print(df_mf_mf.head(1))

cc = df_go_terms.query('go_term_type=="cellular_component"')
df_cc_cc = pd.merge(df_go_edges, cc, 'inner', left_on='x', right_on='go_term_id')
df_cc_cc = df_cc_cc.rename(columns={'go_term_id':'x_id','go_term_name':'x_name','go_term_type':'x_type'})
df_cc_cc = pd.merge(df_cc_cc, cc, 'inner', left_on='y', right_on='go_term_id')
df_cc_cc = df_cc_cc.rename(columns={'go_term_id':'y_id','go_term_name':'y_name','go_term_type':'y_type'})
df_cc_cc['relation'] = 'cellcomp_cellcomp'
df_cc_cc['display_relation'] = 'parent-child'
df_cc_cc['x_source'] = 'GO'
df_cc_cc['y_source'] = 'GO'
df_cc_cc = clean_edges(df_cc_cc)
print(df_cc_cc.head(1))

### Go protein interactions (Gene2GO)

df_prot_path = pd.merge(df_gene2go, df_go_terms, 'inner', 'go_term_id').rename(columns={'go_term_type_x':'go_term_type'})
df_prot_path = pd.merge(df_prot_path, df_prot_names, 'left', left_on='ncbi_gene_id', right_on='ncbi_id')
df_prot_path = df_prot_path.rename(columns={'ncbi_gene_id':'x_id', 'symbol':'x_name', 
                             'go_term_id':'y_id','go_term_name':'y_name', 'go_term_type':'y_type'})
df_prot_path['x_type'] = 'gene/protein'
df_prot_path['x_source'] = 'NCBI'
df_prot_path['y_source'] = 'GO'
df_prot_path = df_prot_path.get(['x_id','x_type', 'x_name', 'x_source','y_id','y_type', 'y_name', 'y_source'])

df_prot_mf = df_prot_path.query('y_type=="molecular_function"').copy()
df_prot_mf['relation'] = 'molfunc_protein'
df_prot_mf['display_relation'] = 'interacts with'
df_prot_mf = clean_edges(df_prot_mf)
print(df_prot_mf.head(1))

df_prot_cc = df_prot_path.query('y_type=="cellular_component"').copy()
df_prot_cc['relation'] = 'cellcomp_protein'
df_prot_cc['display_relation'] = 'interacts with'
df_prot_cc = clean_edges(df_prot_cc)
print(df_prot_cc.head(1))

df_prot_bp = df_prot_path.query('y_type=="biological_process"').copy()
df_prot_bp['relation'] = 'bioprocess_protein'
df_prot_bp['display_relation'] = 'interacts with'
df_prot_bp = clean_edges(df_prot_bp)
print(df_prot_bp.head(1))

## Exposure
### Exposure protein interactions (CTD)

df_exp_prot = df_exposures.get(['exposurestressorname', 'exposurestressorid','exposuremarker', 'exposuremarkerid'])
df_exp_prot = df_exp_prot.loc[df_exp_prot.get(['exposuremarkerid']).dropna().index, :]

gene_row_index = []
for idx, data in df_exp_prot.iterrows():
    if data.exposuremarkerid.isnumeric(): 
        gene_row_index.append(idx)

df_exp_prot = df_exp_prot.loc[gene_row_index, :].astype({'exposuremarkerid': 'int'}).astype({'exposuremarkerid': 'str'})
df_exp_prot = pd.merge(df_exp_prot, df_prot_names, 'left', left_on='exposuremarkerid', right_on='ncbi_id')

df_exp_prot = df_exp_prot.rename(columns={'exposurestressorid':'x_id', 'exposurestressorname':'x_name', 'ncbi_id':'y_id', 'symbol':'y_name'})
df_exp_prot['x_type'] = 'exposure'
df_exp_prot['x_source'] = 'CTD'
df_exp_prot['y_type'] = 'gene/protein'
df_exp_prot['y_source'] = 'NCBI'
df_exp_prot['relation'] = 'exposure_protein'
df_exp_prot['display_relation'] = 'interacts with'
df_exp_prot = clean_edges(df_exp_prot)
print(df_exp_prot.head(1))

### Exposure disease interactions (CTD)

df_exp_dis = df_exposures.get(['exposurestressorname', 'exposurestressorid','diseasename', 'diseaseid'])
df_exp_dis = df_exp_dis.loc[df_exp_dis.get(['diseaseid']).dropna().index, :]
df_exp_dis = pd.merge(df_exp_dis, df_mondo_xref.query('ontology=="MESH"'), 'left', left_on='diseaseid', right_on='ontology_id')
df_exp_dis = pd.merge(df_exp_dis, df_mondo_terms, 'left', left_on='mondo_id', right_on= 'id')

df_exp_dis = df_exp_dis.rename(columns={'exposurestressorid':'x_id', 'exposurestressorname':'x_name', 'mondo_id':'y_id', 'name':'y_name'})
df_exp_dis['x_type'] = 'exposure'
df_exp_dis['x_source'] = 'CTD'
df_exp_dis['y_type'] = 'disease'
df_exp_dis['y_source'] = 'MONDO'
df_exp_dis['relation'] = 'exposure_disease'
df_exp_dis['display_relation'] = 'linked to'
df_exp_dis = clean_edges(df_exp_dis)
print(df_exp_dis.head(1))

### Exposure exposure interactions (CTD)

exposures = np.unique(df_exposures.get('exposurestressorid').values)
df_exp_exp = df_exposures.query('exposuremarkerid in @exposures')

df_exp_exp = df_exp_exp.get(['exposurestressorname', 'exposurestressorid','exposuremarker', 'exposuremarkerid'])
df_exp_exp = df_exp_exp.loc[df_exp_exp.get(['exposuremarkerid']).dropna().index, :]
df_exp_exp = df_exp_exp.drop_duplicates()

df_exp_exp = df_exp_exp.rename(columns={'exposurestressorid':'x_id', 'exposurestressorname':'x_name', 'exposuremarker':'y_name', 'exposuremarkerid':'y_id'})
df_exp_exp['x_type'] = 'exposure'
df_exp_exp['x_source'] = 'CTD'
df_exp_exp['y_type'] = 'exposure'
df_exp_exp['y_source'] = 'CTD'
df_exp_exp['relation'] = 'exposure_exposure'
df_exp_exp['display_relation'] = 'parent-child'
df_exp_exp = clean_edges(df_exp_exp)
print(df_exp_exp.head(1))

### Exposure pathway interactions (CTD)
# phenotypes are actually pathways 

df_exp_path = df_exposures.get(['exposurestressorname', 'exposurestressorid','phenotypename', 'phenotypeid'])
df_exp_path = df_exp_path.loc[df_exp_path.get(['phenotypeid']).dropna().index, :]
df_exp_path.loc[:, 'phenotypeid'] = [str(int(x.split(':')[1])) for x in df_exp_path.get(['phenotypeid']).values.reshape(-1)]
df_exp_path = df_exp_path.drop_duplicates()
df_exp_path = pd.merge(df_exp_path, df_go_terms, 'inner', left_on='phenotypeid', right_on='go_term_id')
df_exp_path = df_exp_path.rename(columns={'exposurestressorid':'x_id', 'exposurestressorname':'x_name', 
                                          'go_term_id':'y_id', 'go_term_name':'y_name', 'go_term_type':'y_type'})
df_exp_path['x_type'] = 'exposure'
df_exp_path['x_source'] = 'CTD'
df_exp_path['y_source'] = 'GO'

df_exp_bp = df_exp_path.query('y_type=="biological_process"').copy()
df_exp_bp['relation'] = 'exposure_bioprocess'
df_exp_bp['display_relation'] = 'interacts with'
df_exp_bp = clean_edges(df_exp_bp)
print(df_exp_bp.head(1))

df_exp_mf = df_exp_path.query('y_type=="molecular_function"').copy()
df_exp_mf['relation'] = 'exposure_molfunc'
df_exp_mf['display_relation'] = 'interacts with'
df_exp_mf = clean_edges(df_exp_mf)
print(df_exp_mf.head(1))

df_exp_cc = df_exp_path.query('y_type=="cellular_component"').copy()
df_exp_cc['relation'] = 'exposure_cellcomp'
df_exp_cc['display_relation'] = 'interacts with'
df_exp_cc = clean_edges(df_exp_cc)
print(df_exp_cc.head(1))

## Anatomy
### Anatomy anatomy interactions (UBERON) 

df_ana_ana = pd.merge(df_uberon_is_a, df_uberon_terms, 'left', left_on='id', right_on='id')
df_ana_ana = df_ana_ana.rename(columns={'id':'x_id', 'name':'x_name'})
df_ana_ana = pd.merge(df_ana_ana, df_uberon_terms, 'left', left_on='is_a', right_on='id')
df_ana_ana = df_ana_ana.rename(columns={'id':'y_id', 'name':'y_name'})
df_ana_ana['x_type'] = 'anatomy'
df_ana_ana['x_source'] = 'UBERON'
df_ana_ana['y_type'] = 'anatomy'
df_ana_ana['y_source'] = 'UBERON'
df_ana_ana['relation'] = 'anatomy_anatomy'
df_ana_ana['display_relation'] = 'parent-child'
df_ana_ana = clean_edges(df_ana_ana)
print(df_ana_ana.head(1))

### Anatomy Protein (BGEE)

df_bgee = pd.merge(df_bgee, df_prot_names, 'inner', left_on='gene_name', right_on='symbol')
df_bgee = df_bgee.rename(columns={'ncbi_id':'x_id', 'symbol':'x_name', 
                                  'anatomy_id':'y_id', 'anatomy_name':'y_name'})
df_bgee['x_source'] = 'NCBI'
df_bgee['x_type'] = 'gene/protein'
df_bgee['y_source'] = 'UBERON'
df_bgee['y_type'] = 'anatomy'

df_ana_prot_pos = df_bgee.query('expression=="present"').copy()
df_ana_prot_pos['relation'] = 'anatomy_protein_present'
df_ana_prot_pos['display_relation'] = 'expression present'
df_ana_prot_pos = clean_edges(df_ana_prot_pos)
print(df_ana_prot_pos.head(1))

df_ana_prot_neg = df_bgee.query('expression=="absent"').copy()
df_ana_prot_neg['relation'] = 'anatomy_protein_absent'
df_ana_prot_neg['display_relation'] = 'expression absent'
df_ana_prot_neg = clean_edges(df_ana_prot_neg)
print(df_ana_prot_neg.head(1))

## Pathways

df_path_path = pd.merge(df_reactome_rels, df_reactome_terms, 'inner', left_on='reactome_id_1', right_on='reactome_id')
df_path_path = df_path_path.rename(columns={'reactome_id': 'x_id', 'reactome_name':'x_name'})
df_path_path = pd.merge(df_path_path, df_reactome_terms, 'inner', left_on='reactome_id_2', right_on='reactome_id')
df_path_path = df_path_path.rename(columns={'reactome_id': 'y_id', 'reactome_name':'y_name'})

df_path_path['x_source'] = 'REACTOME'
df_path_path['x_type'] = 'pathway'
df_path_path['y_source'] = 'REACTOME'
df_path_path['y_type'] = 'pathway'
df_path_path['relation'] = 'pathway_pathway'
df_path_path['display_relation'] = 'parent-child'
df_path_path = clean_edges(df_path_path)
print(df_path_path.head(1))

### Pathway protein interactions

df_path_prot = pd.merge(df_reactome_ncbi, df_prot_names, 'inner', 'ncbi_id')

df_path_prot = df_path_prot.rename(columns={'ncbi_id': 'x_id', 'symbol':'x_name', 
                                            'reactome_id': 'y_id', 'reactome_name':'y_name'})
df_path_prot['x_source'] = 'NCBI'
df_path_prot['x_type'] = 'gene/protein'
df_path_prot['y_source'] = 'REACTOME'
df_path_prot['y_type'] = 'pathway'
df_path_prot['relation'] = 'pathway_protein'
df_path_prot['display_relation'] = 'interacts with'
df_path_prot = clean_edges(df_path_prot)
print(df_path_prot.head(1))

# Compiling knowledge graph

kg = pd.concat([df_prot_drug, df_drug_dis, df_drug_drug, df_prot_phe,
                df_phe_phe, df_dis_phe_neg, df_dis_phe_pos, df_prot_dis, df_dis_dis, 
                df_drug_effect, df_bp_bp, df_mf_mf, df_cc_cc, df_prot_mf, 
                df_prot_cc, df_prot_bp, df_exp_prot, df_exp_dis, df_exp_exp, 
                df_exp_bp, df_exp_mf, df_exp_cc, df_path_path, df_path_prot,
                df_ana_ana, df_ana_prot_pos, df_ana_prot_neg]) #28
kg = kg.drop_duplicates()
kg_rev = kg.copy().rename(columns={'x_id':'y_id','x_type':'y_type', 'x_name':'y_name', 'x_source':'y_source',
                            'y_id':'x_id','y_type':'x_type', 'y_name':'x_name', 'y_source':'x_source' }) #add reverse edges
kg = pd.concat([kg, kg_rev])
kg = kg.drop_duplicates()
kg = kg.dropna()
# remove self loops from edges 
kg = kg.query('not ((x_id == y_id) and (x_type == y_type) and (x_source == y_source) and (x_name == y_name))')
print(kg.head())

kg.to_csv(SAVE_PATH+'kg_raw.csv', index=False)

# Get giant component

nodes = pd.concat([kg.get(['x_id','x_type', 'x_name','x_source']).rename(columns={'x_id':'node_id', 'x_type':'node_type', 'x_name':'node_name','x_source':'node_source'}), 
                   kg.get(['y_id','y_type', 'y_name','y_source']).rename(columns={'y_id':'node_id', 'y_type':'node_type', 'y_name':'node_name','y_source':'node_source'})])
nodes = nodes.drop_duplicates().reset_index().drop('index',axis=1).reset_index().rename(columns={'index':'node_idx'})

edges = pd.merge(kg, nodes, 'left', left_on=['x_id','x_type', 'x_name','x_source'], right_on=['node_id','node_type','node_name','node_source'])
edges = edges.rename(columns={'node_idx':'x_idx'})
edges = pd.merge(edges, nodes, 'left', left_on=['y_id','y_type', 'y_name','y_source'], right_on=['node_id','node_type','node_name','node_source'])
edges = edges.rename(columns={'node_idx':'y_idx'})
edges = edges.get(['relation', 'display_relation','x_idx', 'y_idx'])
edges['combine_idx'] = edges['x_idx'].astype(str) + '-' + edges['y_idx'].astype(str)

edge_index = edges.get(['x_idx', 'y_idx']).values.T

graph = ig.Graph()
graph.add_vertices(list(range(nodes.shape[0])))
graph.add_edges([tuple(x) for x in edge_index.T])

graph = graph.as_undirected(mode='collapse')

c = graph.components(mode='strong')
giant = c.giant()

print('Nodes: %d' % giant.vcount())
print('Edges: %d' % giant.ecount())

assert not giant.is_directed()
assert giant.is_connected()

giant_nodes = giant.vs['name']
new_nodes = nodes.query('node_idx in @giant_nodes')
assert new_nodes.shape[0] == giant.vcount()

new_edges = edges.query('x_idx in @giant_nodes and y_idx in @giant_nodes').copy()
assert new_edges.shape[0] == giant.ecount()

new_kg = pd.merge(new_edges, new_nodes, 'left', left_on='x_idx', right_on='node_idx')
new_kg = new_kg.rename(columns={'node_id':'x_id', 'node_type':'x_type', 'node_name':'x_name','node_source':'x_source'}) 
new_kg = pd.merge(new_kg, new_nodes, 'left', left_on='y_idx', right_on='node_idx')
new_kg = new_kg.rename(columns={'node_id':'y_id', 'node_type':'y_type', 'node_name':'y_name','node_source':'y_source'}) 
new_kg = clean_edges(new_kg)

kg = new_kg.copy()
kg.to_csv(SAVE_PATH+'kg_giant.csv', index=False)