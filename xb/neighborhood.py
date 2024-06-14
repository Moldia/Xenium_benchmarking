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


def nhood_squidpy(adata,sample_key='sample',radius=50,cluster_key='leiden',save=True,plot_path='./',cmap='inferno',vmax=None,vmin=None):
    """ Compute neighborhood enrichment based on Squidpy's function

    Parameters:
    adata (AnnData): AnnData object with the cells of the experiment. 
    sample_key (str): name of the column where the sample each cell belongs to is specify. It should be a column present in adata.obs
    radius (int): radius to consider when compuing the spatial neighbors, specified in the scale that adata.obsm['spatial'] is in (typically um)
    cluster_key (str): name of the column where the cell type of each cell is specified. The neighborhood enrichment will be computed based on this groups 
    save (Boolean): specify whether the resulting plot should be saved in the paths specified in 'plot_path' or not
    cmap (str): name of the colormap used to plot the neighborhood enrichment plot
    vmax (int): maximum value to show in the neighborhood enrcihment plot 
    vmin (int): minimum value to show in the neighborhood enrichment plot

    Returns:
    adata1: AnnData object with the neighborhood enrichment scores computed

    """
    
    anndata_list=[]
    plt.style.use('default')
    for sample in adata.obs[sample_key].unique():
        adata_copy_int = adata[adata.obs[sample_key]==sample]
        sq.gr.spatial_neighbors(adata_copy_int,radius=radius,coord_type ='generic')
        anndata_list.append(adata_copy_int) 
    adata1=sc.concat(anndata_list,join='outer',pairwise=True) 
    sq.gr.nhood_enrichment(adata1, cluster_key=cluster_key)
    sq.pl.nhood_enrichment(adata1, cluster_key=cluster_key, method="single", cmap=cmap, vmin=vmin, vmax=vmax,show=False)
    if save==True:
        plt.savefig(plot_path+'neighborhood_enrichment_'+str(cluster_key)+'_'+str(radius)+'.pdf')
    return adata1



