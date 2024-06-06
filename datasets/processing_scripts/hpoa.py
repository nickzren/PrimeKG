import os
import numpy as np
import pandas as pd

# Reading annotations file 
with open('data/input/hpo/phenotype.hpoa', 'r') as f: 
    data = f.readlines()

# Removing the first four lines
data_str = ''.join(data[4:])

# Writing the formatted data to a new file
with open('data/input/hpo/phenotype_formatted.hpoa', 'w') as f:
    f.write(data_str)

# Specify the data types for the problematic columns
dtype_spec = {
    2: 'str',
    6: 'str',
    8: 'str',
    9: 'str'
}

# Reading the CSV file with the specified dtype
df = pd.read_csv('data/input/hpo/phenotype_formatted.hpoa', sep='\t', dtype=dtype_spec)

# Processing the data
for x in df.itertuples():
    ont, ont_id = x.database_id.split(':')
    df.loc[x.Index, 'disease_ontology'] = ont
    df.loc[x.Index, 'disease_ontology_id'] = ont_id
    df.loc[x.Index, 'hp_id'] = str(int(x.hpo_id.split(':')[1]))

# Writing the negative and positive phenotype data to CSV files
df.query('qualifier=="NOT"').get(['hp_id', 'disease_ontology', 'disease_ontology_id']).drop_duplicates().to_csv('data/output/hpo/disease_phenotype_neg.csv', index=False)
df.query('qualifier!="NOT"').get(['hp_id', 'disease_ontology', 'disease_ontology_id']).drop_duplicates().to_csv('data/output/hpo/disease_phenotype_pos.csv', index=False)
