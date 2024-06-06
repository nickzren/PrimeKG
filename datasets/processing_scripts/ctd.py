import numpy as np
import pandas as pd

# Read the lines from the file
with open('data/input/ctd/CTD_exposure_events.csv', 'r') as f:
    lines = f.readlines()

# Extract the column names
field_next = False
for line in lines:
    if line.startswith('# Fields'):
        field_next = True
        continue
    if field_next:
        fields = line
        break
cols = fields[2:-2].split(',')

# Write the column names and data to a new CSV file in the input directory
with open('data/input/ctd/exposure_data.csv', 'w') as f:
    f.write(','.join(cols) + '\n')
    for line in lines:
        if not line.startswith('#'):
            f.write(line + '\n')

# Specify the data types for columns, setting column 16 (17th column) to 'str'
dtype_spec = {col: 'str' for col in range(len(cols))}

# Read the CSV file with the specified dtype from the input directory
df = pd.read_csv('data/input/ctd/exposure_data.csv', index_col=False, dtype=dtype_spec, low_memory=False)

# Save the processed data to the output directory
df.to_csv('data/output/ctd/exposure_data.csv', index=False)