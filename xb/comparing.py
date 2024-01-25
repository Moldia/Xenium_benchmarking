import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
import shutil
import os.path
from scipy.io import mmread
import tifffile as tf
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from anndata import AnnData
import json
from sklearn.metrics import mutual_info_score
from sklearn.metrics import fowlkes_mallows_score

def combine_med(medians,tag):
    mm=pd.DataFrame(medians,columns=['ratio'])
    mm['method']=tag
    return mm
def median_calculator(adata_dict,df_filt):
    gg=['anno_scRNAseq','anno_CosMx', 'anno_Hybriss', 'anno_MERFISH', 'anno_ResolveBio', 'anno_Vizgen', 'anno_Xenium']
    genes_s  = {
        'anno_scRNAseq': [],
         'anno_starmap': [],
         'anno_allen_smfish': [],
         'anno_MERFISH': [],
        'anno_Hybriss': [],
         'anno_osmfish':[],
         'anno_seqFISH': [],
         'anno_exseq': [],
         'anno_Vizgen':[],
         'anno_baristaseq': [],
         'anno_Xenium': [],
        'anno_ResolveBio': [],
        'anno_CosMx': [],
    }
    means  = {
        'anno_scRNAseq': [],
         'anno_starmap': [],
         'anno_allen_smfish': [],
         'anno_MERFISH': [],
        'anno_Hybriss': [],
         'anno_osmfish':[],
         'anno_seqFISH': [],
         'anno_exseq': [],
         'anno_Vizgen':[],
         'anno_baristaseq': [],
         'anno_Xenium': [],
        'anno_ResolveBio': [],
        'anno_CosMx': [],
    }
    for i,gene in tqdm(enumerate(df_filt.index)):
        RNA_S=[]
        print(gene)
        if np.sum(adata_dict['anno_scRNAseq'][:, adata_dict['anno_scRNAseq'].var.index == (gene)].X>minreads)>0:
            for key in gg:
                if key == 'anno_scRNAseq':
                    rna_subset = adata_dict[key][:, adata_dict[key].var.index == (gene)]
                    print(rna_subset.shape)
                    rna_median = np.median(list(rna_subset.X[rna_subset.X > minreads]))
                else:
                    if gene in adata_dict[key].var.index:
                        subset = adata_dict[key][:, adata_dict[key].var.index == (gene)]
                        median_exp = np.median(subset.X[subset.X > minreads])
                        #median_exp/rna_median
                        df_int = pd.DataFrame(subset.X[subset.X > minreads])
                        df_int['method'] = key
                        df_int[0].mean()/rna_median
                        RNA_S.append(df_int)
                        means[key].append(df_int[0].mean()/rna_median)
                        genes_s[key].append(gene)
                    if len(RNA_S) < 1:
                        continue
                    else:
                        RNA_df = pd.concat(RNA_S)
                        RNA_df=RNA_df.reset_index()

                    #means[key].append(df_int[0].mean()/rna_median)

            RNA_df['color'] = RNA_df.method.map(color_dicitonary)
            #UNCOMMENT TO GET THE PLOT
        #    sns.displot(data=RNA_df,x=RNA_df[0],hue='method',kind='ecdf',log_scale=True, palette=list(RNA_df['color'].unique()))
        #    plt.title(gene)
        #    plt.show()
        #UNTIL HERE
    return means,genes_s