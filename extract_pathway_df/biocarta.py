import pandas as pd
import sys
from bs4 import BeautifulSoup
import requests

def getDescription(pathway_id):
    # example url: http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=BIOCARTA_CK1_PATHWAY
    base_url = "http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=" + pathway_id
    response = requests.get(base_url)
    soup = BeautifulSoup(response.text, 'html.parser')
    description_td = [t.parent.findNext('td') for t in soup.findAll(text='Brief description')][0]
    return description_td.text
    
getDescription('BIOCARTA_CK1_PATHWAY')

# path for biocarta.gmt
panter_path = sys.argv[1]
pathways = []
content = None

with open(panter_path) as f:
    content = f.read()

for line in content.split('\n'):
    parts = line.split('\t')
    pathway_id = parts[0].split('%')[-1]
    try:
        pathway_name = getDescription(pathway_id)
    except:
        pathway_name = pathway_id
        
    # PATHWAY_NAME%PANTHER PATHWAY%P04373   NAME    ENT1    ENT2 ...
    for entrez_id in parts[2:]:
        if entrez_id:
            pathway_info = {'database': 'biocarta', 'pathway_id': pathway_id, 'pathway_name': pathway_name, 'entrez': entrez_id}
            pathways.append(pathway_info)

df = pd.DataFrame(pathways)
df.to_pickle('biocarta.pkl')



