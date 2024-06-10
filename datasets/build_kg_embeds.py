import os
import torch
from transformers import AutoTokenizer, AutoModel, pipeline
import pandas as pd
import numpy as np
from tqdm import tqdm
import shutil
from scipy.sparse import lil_matrix, save_npz
from sklearn.metrics.pairwise import cosine_similarity
import json

# Constants
INPUT = 'data/input/'
OUTPUT = 'data/output/'
SAVE_PATH = OUTPUT + 'kg/'

kg = pd.read_csv(SAVE_PATH + 'kg_giant.csv', low_memory=False)
disease_nodes = pd.concat([kg.get(['x_id', 'x_type', 'x_name', 'x_source']).rename(columns={'x_id': 'node_id', 'x_type': 'node_type', 'x_name': 'node_name', 'x_source': 'node_source'}),
                           kg.get(['y_id', 'y_type', 'y_name', 'y_source']).rename(columns={'y_id': 'node_id', 'y_type': 'node_type', 'y_name': 'node_name', 'y_source': 'node_source'})])
disease_nodes = disease_nodes.query('node_type=="disease"')
disease_nodes = disease_nodes.drop_duplicates().reset_index().drop('index', axis=1).reset_index().rename(columns={'index': 'node_idx'})

idx2group = {}
with open(SAVE_PATH + 'idx2group.txt', 'r') as f:
    idx2group = json.load(f)

idx2group = {int(k): v for k, v in idx2group.items()}

disease_nodes.loc[:, 'group_name'] = ''

for data in disease_nodes.itertuples():
    if data.node_idx in idx2group.keys():
        disease_nodes.loc[data.Index, 'group_name'] = str(idx2group[data.node_idx])
    else:
        disease_nodes.loc[data.Index, 'group_name'] = data.node_name

disease_group_1 = disease_nodes.get(['group_name']).drop_duplicates().reset_index().rename(columns={'index': 'group_idx'})
disease_group_1.to_csv(SAVE_PATH + 'disease_group_1.csv', index=False)

disease_nodes = pd.merge(disease_nodes, disease_group_1, 'left', 'group_name')
disease_nodes.to_csv(SAVE_PATH + 'disease_nodes.csv', index=False)


input_text = list(disease_group_1.get('group_name').values)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
model_name = 'emilyalsentzer/Bio_ClinicalBERT'
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)
model = model.to(device)
model.eval()

print("Model set to evaluation mode")

def batch(iterable, batch_size=4, return_idx=True):
    l = len(iterable)
    for ndx in range(0, l, batch_size):
        if return_idx:
            yield (ndx, min(ndx + batch_size, l))
        else:
            yield iterable[ndx:min(ndx + batch_size, l)]

tmp_dir = 'tmp/'
if os.path.isdir(tmp_dir):
    shutil.rmtree(tmp_dir)
os.mkdir(tmp_dir)

batch_size = 32
input_tokens = tokenizer(input_text, padding=True, return_tensors='pt', truncation=True, max_length=512)
total_batches = len(input_text) // batch_size + (1 if len(input_text) % batch_size != 0 else 0)

for i, (start, end) in tqdm(enumerate(batch(input_text, batch_size)), total=total_batches):
    input_ids = input_tokens['input_ids'][start:end, :].to(device)
    attention_mask = input_tokens['attention_mask'][start:end, :].to(device)
    with torch.no_grad():
        outputs = model(input_ids=input_ids, attention_mask=attention_mask)
        embeds = torch.mean(outputs[0], dim=1)
    np.save(tmp_dir + str(i) + '.npy', embeds.cpu().numpy())

embeds = []
for i, _ in tqdm(enumerate(batch(input_text, batch_size)), total=total_batches):
    x = np.load(tmp_dir + str(i) + '.npy')
    embeds.append(x)
embeds = np.concatenate(embeds)

np.save(SAVE_PATH + 'kg_disease_bert_embeds.npy', embeds)
if os.path.isdir(tmp_dir):
    shutil.rmtree(tmp_dir)
