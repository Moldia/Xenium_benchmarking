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
import squidpy as sq


def nhood_squidpy(adata1,sample_key='sample',radius=50.0,cluster_key='leiden',save=True,plot_path='./'):
    anndata_list=[]
    plt.style.use('default')
    for sample in adata.obs[sample_key].unique():
        adata_copy_int = adata[adata.obs[sample_key]==sample]
        sq.gr.spatial_neighbors(adata_copy_int,radius=radius,coord_type ='generic')
        anndata_list.append(adata_copy_int) 
    adata1=sc.concat(anndata_list,join='outer',pairwise=True) 
    sq.gr.nhood_enrichment(adata1, cluster_key=cluster_key)
    sq.pl.nhood_enrichment(adata1, cluster_key=cluster_key, method="single", cmap="inferno", vmin=-50, vmax=100,show=False)
    if save==True:
        plt.savefig(plot_path+'neighborhood_enrichment_'+str(cluster_key)+'_'+str(radius)+'.pdf')
    return adata1



