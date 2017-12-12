import pandas as pd
import os
import xml.etree.ElementTree
import sys

# Open a file
path = sys.argv[1]
dirs = os.listdir(path)

pathways = []

for file_path in dirs:
    if '.gpml' in file_path:
        root = xml.etree.ElementTree.parse(file_path).getroot()
        

        path_way_name = root.get('Name')
        file_name = file_path.replace('.gpml', '')
        pathway_id = file_name.split('_')[-1]
        

        for datanode in root.findall('{http://pathvisio.org/GPML/2013a}DataNode'):
            atype = datanode.find('{http://pathvisio.org/GPML/2013a}Xref')
            if 'Entrez Gene' in atype.get('Database'):
                gene = atype.get('ID')
                if gene:
                    pathway_info  = {'database': 'wikipathway', 'pathway_id': pathway_id, 'pathway_name': path_way_name, 'entrez': gene}
                    pathways.append(pathway_info)
        
df = pd.DataFrame.from_dict(pathways)
df.to_pickle('wikipathway.pkl')
