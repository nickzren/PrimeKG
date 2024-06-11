import pandas as pd
import re
import json

# Constants
INPUT = 'data/input/'
OUTPUT = 'data/output/'
SAVE_PATH = OUTPUT + 'kg/'

kg = pd.read_csv(SAVE_PATH+'kg_giant.csv', low_memory=False)

disease_nodes = pd.concat([kg.get(['x_id','x_type', 'x_name','x_source']).rename(columns={'x_id':'node_id', 'x_type':'node_type', 'x_name':'node_name','x_source':'node_source'}), 
                   kg.get(['y_id','y_type', 'y_name','y_source']).rename(columns={'y_id':'node_id', 'y_type':'node_type', 'y_name':'node_name','y_source':'node_source'})])
disease_nodes = disease_nodes.query('node_type=="disease"')
disease_nodes = disease_nodes.drop_duplicates().reset_index().drop('index',axis=1).reset_index().rename(columns={'index':'node_idx'})

seen = set()
idx2group = {}
no = set()

def isroman(s):
    return bool(re.search(r"^M{0,3}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$",s))

def issingleletter(s): 
    if len(s)>1: return False

def same_words(s1, s2): 
    for word in s1.lower().split(' '): 
        word = word.split(',')[0]
        if word!='type' and word!='(disease)' and word not in s2.lower(): 
            return False 
    for word in s2.lower().split(' '): 
        word = word.split(',')[0]
        if word!='type' and word!='(disease)' and word not in s1.lower(): 
            return False
    return True

for i in range(disease_nodes.shape[0]):
    i_name = disease_nodes.loc[i, 'node_name']
    i_idx = disease_nodes.loc[i, 'node_idx']
    for w in ['monosomy','disomy', 'trisomy', 'trisomy/tetrasomy', 'chromosome']: 
        if w in i_name: 
            no.add(i_idx)

for i in range(disease_nodes.shape[0]):
    i_idx = disease_nodes.loc[i, 'node_idx']
    if i_idx in seen: continue 
    if i_idx in no: continue 
    i_name = disease_nodes.loc[i, 'node_name']
    i_split = i_name.split(' ')
    end = i_split[-1]
    if len(end)<=2 or end.isnumeric() or isroman(end):  
        main_text = ' '.join(i_split[:-1])
        matches = [i_name]
        matches_idx = [i_idx]
        match_found = False
        numeric = True
        for j in range(disease_nodes.shape[0]):
            j_idx = disease_nodes.loc[j, 'node_idx']
            j_name = disease_nodes.loc[j, 'node_name']
            m = ' '.join(j_name.split(' ')[:-1])
            if m.lower() == main_text.lower() or same_words(m, main_text): 
                matches.append(j_name)
                matches_idx.append(j_idx)
                match_found = True
        if match_found:
            matches_idx = list(set(matches_idx))
            matches = list(set(matches))
            if len(matches) <= 1: continue 
            if main_text.endswith('type'): 
                main_text = main_text[:-4]
            if main_text.endswith(','): 
                main_text = main_text[:-1]
            if main_text.endswith(' '): 
                main_text = main_text[:-1]
            for x in matches_idx: 
                seen.add(x)
                idx2group[int(x)] = main_text

with open(SAVE_PATH + 'idx2group.txt', 'w') as f:
    json.dump(idx2group, f, ensure_ascii=False, indent=4)