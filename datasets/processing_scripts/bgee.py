import numpy as np
import pandas as pd

# Specify the data type for column 22 (assuming it's the 23rd column since indexing starts from 0)
dtype_spec = {22: 'str'} 

# Read the CSV file with the specified dtype
df = pd.read_csv('data/input/bgee/Homo_sapiens_expr_advanced.tsv', sep='\t', dtype=dtype_spec)

# Filter and rename columns
df = df[df.get('Anatomical entity ID').str.startswith('UBERON')]
df = df.get(['Gene ID', 'Gene name', 'Anatomical entity ID', 'Anatomical entity name', 'Expression', 'Call quality', 'Expression rank'])

df = df.rename(columns={'Gene ID':'gene_id',
                        'Gene name':'gene_name',
                        'Anatomical entity ID':'anatomy_id',
                        'Anatomical entity name':'anatomy_name',
                        'Expression':'expression',
                        'Call quality':'call_quality',
                        'Expression rank': 'expression_rank'})

# Apply filters
df = df.query('call_quality=="gold quality"')
df = df.query('expression_rank<25000')

# Updated PrimeKG: remove rows that specify cell type within a particular tissue
df = df[~df['anatomy_id'].str.contains("∩")]

# Process anatomy_id
df.loc[:, 'anatomy_id'] = [str(int(x.split(':')[1])) for x in df.get(['anatomy_id']).values.reshape(-1)]

# Save the processed data to CSV
df.to_csv('data/output/bgee/anatomy_gene.csv', index=False)
