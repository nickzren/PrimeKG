import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity

# Constants
INPUT = 'data/input/'
OUTPUT = 'data/output/'
SAVE_PATH = OUTPUT + 'kg/'

kg = pd.read_csv(SAVE_PATH + 'kg_giant.csv', low_memory=False)
grouped_diseases = pd.read_csv(SAVE_PATH+'kg_grouped_diseases.csv').astype({'node_id':str})
group_col = 'group_name_bert'

groups = grouped_diseases.groupby(group_col).count().query('node_id>1').index.values
set_groups = set(groups)

id_col = group_col.replace('name','id')
group_map = pd.DataFrame(columns=[id_col, group_col])
group_map.loc[:, group_col] = groups

grouped_diseases = grouped_diseases.query('{} in @set_groups'.format(group_col))

for g, data in grouped_diseases.groupby(group_col): 
    if g in set_groups:
        x = '_'.join(list(data.get('node_id').values))
        i = group_map.query('{}==@g'.format(group_col)).index[0]
        group_map.loc[i, id_col] = x
        
grouped_diseases = pd.merge(grouped_diseases, group_map)
grouped_diseases.to_csv(SAVE_PATH+'kg_grouped_diseases_bert_map.csv', index=False)

kg_x_dis = kg.query('x_type=="disease" and x_source=="MONDO"')
kg_y_dis = kg.query('y_type=="disease" and y_source=="MONDO"')

for idx, data in tqdm(grouped_diseases.iterrows(), total=grouped_diseases.shape[0]): 
    x_index = kg_x_dis.query('x_id==@data.node_id and x_name==@data.node_name').index.values
    kg.loc[x_index, 'x_id'] = data.get(id_col)
    kg.loc[x_index, 'x_name'] = data.get(group_col)
    kg.loc[x_index, 'x_source'] = 'MONDO_grouped'

    y_index = kg_y_dis.query('y_id==@data.node_id and y_name==@data.node_name').index.values
    kg.loc[y_index, 'y_id'] = data.get(id_col)
    kg.loc[y_index, 'y_name'] = data.get(group_col)
    kg.loc[y_index, 'y_source'] = 'MONDO_grouped'

kg = kg.drop_duplicates()
kg_rev = kg.copy().rename(columns={'x_id':'y_id','x_type':'y_type', 'x_name':'y_name', 'x_source':'y_source',
                                   'y_id':'x_id','y_type':'x_type', 'y_name':'x_name', 'y_source':'x_source' }) #add reverse edges
kg = pd.concat([kg, kg_rev])
kg = kg.drop_duplicates()
kg = kg.dropna()
# remove self loops from edges 
kg = kg.query('not ((x_id == y_id) and (x_type == y_type) and (x_source == y_source) and (x_name == y_name))')

kg.to_csv(SAVE_PATH+'kg_grouped.csv', index=False)

kg = pd.read_csv(SAVE_PATH+'kg_grouped.csv', low_memory=False)

# nodes file 
nodes = pd.concat([kg.get(['x_id','x_type', 'x_name','x_source']).rename(columns={'x_id':'node_id', 'x_type':'node_type', 'x_name':'node_name', 'x_source':'node_source'}), 
                   kg.get(['y_id','y_type', 'y_name','y_source']).rename(columns={'y_id':'node_id', 'y_type':'node_type', 'y_name':'node_name', 'y_source':'node_source'})])
nodes = nodes.drop_duplicates().reset_index().drop('index',axis=1).reset_index().rename(columns={'index':'node_index'})

# assign index 
kg = pd.merge(kg, nodes.rename(columns={'node_index':'x_index',
                                        'node_id':'x_id',
                                        'node_type':'x_type',
                                        'node_name':'x_name',
                                        'node_source':'x_source'}), 'left').dropna()
kg = pd.merge(kg, nodes.rename(columns={'node_index':'y_index',
                                        'node_id':'y_id',
                                        'node_type':'y_type',
                                        'node_name':'y_name',
                                        'node_source':'y_source'}), 'left').dropna()
kg = kg.get(['relation', 'display_relation', 'x_index', 'x_id', 'x_type', 'x_name', 'x_source',
       'y_index', 'y_id', 'y_type', 'y_name', 'y_source'])

# edges file 
edges = kg.get(['relation', 'display_relation', 'x_index', 'y_index']).copy()

kg.to_csv(SAVE_PATH+'kg.csv', index=False)
nodes.to_csv(SAVE_PATH+'nodes.csv', index=False)
edges.to_csv(SAVE_PATH+'edges.csv', index=False)

def kg_describe(df, by, count_col): 
    df = df.groupby(by).count().sort_values(by=count_col, ascending=False).rename(columns={count_col: 'count'}).get(['count'])
    total = np.sum(df.get('count').values)
    df = df.eval('percent = 100*count/@total')
    total_row = pd.DataFrame(df.sum(0)).T
    total_row.index = ['total']
    df = pd.concat([df, total_row])
    df['count'] = df.get(['count']).astype('int')
    df['percent'] = df.get(['percent']).round(1)
    print(df)
    return df

# Call the function and print the descriptions
print("Node Type Description:")
kg_describe(nodes, 'node_type', 'node_index')

print("\nRelation Description:")
kg_describe(edges, 'relation', 'x_index')