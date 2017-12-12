import pandas as pd
from kegg.lib import get_kegg_data
import numpy as np

kegg_df = get_kegg_data("hsa", True)
kegg_df.columns = ['pathway_id', 'entrez', 'pathway_name']
kegg_df['database'] = pd.Series(np.repeat('kegg', kegg_df.shape[0]), index=kegg_df.index)
kegg_df.to_pickle('kegg.pkl')


