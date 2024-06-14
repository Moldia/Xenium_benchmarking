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
from tqdm import tqdm

def combine_med(medians,tag):
    """ Combine precomputed medians into a single dataframe
   
    Parameters:
    medians(list): list of precomputed medians of expression
    tag(str): tag to be added as a column to the list of medians. In here, this is the methods the medians where computed from
    Returns:
    mm(DataFrame): formated medians into a DataFrame

   """
    mm=pd.DataFrame(medians,columns=['ratio'])
    mm['method']=tag
    return mm


def median_calculator(adata_dict,df_filt):
    """ Calculate medians expression for cells profiled with each technology compared to a reference single cell RNAseq dataset
   
    Parameters:
    adata_dict (dict): dictionary including the names of the datasets analyzed as .keys() and AnnData's of each technologies as .values(). It includes a reference scRNAseq dataset in 'anno_scRNAseq'
    df_filt(DataFrame): dataframe including the list of genes to be compared in .index.  
    Returns:
    means(dict): dictionary of means computed with names of the datasets in .keys() and a list of medians computed as .values()
    genes_s(dict):dictionary of gene name of the means computed with names of the datasets in .keys() and a list of neme of the genes that have been used to compute medians computed as .values()

   """
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

            RNA_df['color'] = RNA_df.method.map(color_dicitonary)

    return means,genes_s
