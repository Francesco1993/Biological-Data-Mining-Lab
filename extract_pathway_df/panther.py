import pandas as pd
import sys

panter_path = sys.argv[1]
pathways = []
content = None

with open(panter_path) as f:
    content = f.read()

for line in content.split('\n'):
    parts = line.split('\t')
    pathway_id = parts[0].split('%')[-1]
    pathway_name = parts[0].split('%')[0].lower()
    # PATHWAY_NAME%PANTHER PATHWAY%P04373   NAME    ENT1    ENT2 ...
    for entrez_id in parts[2:]:
        pathway_info = {'database': 'panther', 'pathway_id': pathway_id, 'pathway_name': pathway_name, 'entrez': entrez_id}
        pathways.append(pathway_info)

df = pd.DataFrame.from_dict(pathways)
df.to_pickle('panther.pkl')



