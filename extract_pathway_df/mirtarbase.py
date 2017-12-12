import pandas as pd

df = pd.read_csv("./data/hsa_mirTARbase.csv", sep=",")

a = list(df["miRNA"])
b = list(df["Target Gene (Entrez Gene ID)"])

pathways = []

#PYTHON
support_dictio = {}
for el1, el2 in zip(a,b):
    if not el1 in support_dictio.keys():
        support_dictio[el1] = []
        support_dictio[el1].append(str(el2))
    elif el1 in support_dictio.keys() and str(el2) not in support_dictio[el1]:
        support_dictio[el1].append(str(el2))
    else:
        continue

for key, entrez_list in support_dictio.items():
    for ent in entrez_list:
        pathway_info =  {'database': 'MiRTarBase', 'pathway_id': key, 'pathway_name': key, 'entrez': ent}
        pathways.append(pathway_info)

df = pd.DataFrame(pathways)
df.to_pickle('mirtarbase.pkl')
