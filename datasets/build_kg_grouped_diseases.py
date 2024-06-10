import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity

# Constants
INPUT = 'data/input/'
OUTPUT = 'data/output/'
SAVE_PATH = OUTPUT + 'kg/'

disease_group_1 = pd.read_csv(SAVE_PATH + 'disease_group_1.csv')
disease_nodes = pd.read_csv(SAVE_PATH + 'disease_nodes.csv')

embeds = np.load(SAVE_PATH+'kg_disease_bert_embeds.npy')
cos_sim = cosine_similarity(embeds, embeds)

seen = set()
groups = []
idx2group = {}
no = set()

for i in range(disease_group_1.shape[0]):
    i_name = disease_group_1.loc[i, 'group_name']
    i_idx = disease_group_1.loc[i, 'group_idx']
    for w in ['cardiomyopathy', 'syndrome', 'combined', 'complement', 'deficiency', 
              'factor', 'immunodeficiency', 'monosomy','disomy', 'trisomy', 
              'trisomy/tetrasomy', 'chromosome', 'neuroendocrine tumor', 
              'neuroendocrine neoplasm', 'cancer', 'tumor', 'neoplasm','carcinoma',
              'lymphoma', 'lipoma']: 
        if w in i_name: 
            no.add(i_idx)
            continue
    for w in ['CDG']: 
        if i_name.endswith(w): 
            no.add(i_idx)
            continue
    for w in ['neurodevelopmental disorder', 'glycogen storage disease', 
              'congenital disorder of glycosylation', 'qualitative or quantitative defects']: 
        if i_name.startswith(w): 
            no.add(i_idx)
            continue
            
cutoff = 0.98
for i in range(disease_group_1.shape[0]):
    i_name = disease_group_1.loc[i, 'group_name']
    i_idx = disease_group_1.loc[i, 'group_idx']
    if i_idx in no or i_idx in seen: continue
    x = disease_group_1[cos_sim[i]>cutoff]
    if x.shape[0]>1: 
        main_text = i_name  # Automatically use the first disease group name
        for v in x.get('group_name').values: 
            print(v)
        print(f'Assigned group name: {main_text}')  # Inform which group name was assigned
        for v in x.get('group_idx').values: 
            seen.add(v)
            idx2group[v] = main_text
        g = list(x.get('group_idx').values.reshape(-1))
        groups.append((main_text, g))  # main_text contains group name
    else: 
        no.add(i_idx)
        print('Not added')


disease_group_1.loc[:, 'group_name_2'] = ''
for data in disease_group_1.itertuples(): 
    if data.group_idx in idx2group.keys():
        disease_group_1.loc[data.Index, 'group_name_2'] = idx2group[data.group_idx]
    else: 
        disease_group_1.loc[data.Index, 'group_name_2'] = data.group_name
        
disease_group_2 = disease_group_1.get(['group_name_2']).drop_duplicates().reset_index().rename(columns={'index':'group_idx_2'})

df_disease_group = pd.merge(disease_nodes, disease_group_1, 'left', 'group_name')
df_disease_group = df_disease_group.get(['node_id', 'node_type', 'node_name', 'node_source',
       'group_name', 'group_name_2'])
df_disease_group = df_disease_group.rename(columns={'group_name':'group_name_auto',
        'group_name_2':'group_name_bert'}).astype({'node_id':str})
df_disease_group.to_csv(SAVE_PATH+'kg_grouped_diseases.csv')