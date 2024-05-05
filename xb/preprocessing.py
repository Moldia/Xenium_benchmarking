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
import scipy as sp
from xb.plotting import plot_cell_counts

def main_preprocessing(adata,target_sum=100,mincounts=10,mingenes=3,neigh=15,npc=0,nuc=1,scale=False,hvg=False,default=False,total_clusters=30,norm=True,lg=True):
    """ Preprocess and cluster the cells in adata, given the parameters specified. This function is mainly used for simulating the performance of different preprocessing strategies
   
    Parameters:
    adata (AnnData): AnnData object with the cells of the experiment. 
    norm(boolean): Whether to normalize based cells or not
    target_sum(int or None): Target sum to use if the normalization is done based on library size. None is used for automatic calculation of library size
    lg(boolean): Whether to log-transforms cells.
    mincounts (int): Minimum amount of counts detected in a cell to pass the quality filters
    mingenes (int): Minimum amount of genes expressed in a cell to pass the quality filters
    neigh(int): number of neighbors to used when calculating the nearest neighbors by sc.pp.neighbors()
    npc(int): number of principal components to used when calculating the nearest neighbors by sc.pp.neighbors()
    scale(boolean): whether to scale the data or not
    hvg(boolean): whether to select highly variable genes for further processing or not
    total_clusters (int): number of clusters to obtain in the process of clustering (+-2)
    default(boolean): whether the run is the original one or not
     nuc(int): DEPRECATED. NOT USED IN THIS FUNCTION
    
    
    

    Returns:
    adata: AnnData object with the preprocessed and clustered cells according to the parameters specified

   """
    adata.layers['raw']=adata.X.copy()
    sc.pp.filter_cells(adata,min_counts=mincounts)
    sc.pp.filter_cells(adata,min_genes=mingenes)
    adata.raw=adata
    adata.layers['raw']=adata.X.copy()
    if norm==True:
        sc.pp.normalize_total(adata, target_sum=target_sum)
    if lg==True:
        sc.pp.log1p(adata)
    if hvg==True:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
    if scale==True:
        sc.pp.scale(adata)
    sc.pp.pca(adata)
    print(adata.shape)
    sc.pp.neighbors(adata, n_neighbors=neigh, n_pcs=npc)
    resol=1.4
    sc.tl.leiden(adata,resolution=resol,key_added='leiden_1_4')
    numclust=int(np.max(adata.obs['leiden_1_4'].astype(int)))
    print(numclust)
    targetnum=total_clusters
    if default==True:
        targetnum=numclust
    if (numclust-targetnum)>3:
        i=0
        while (int(numclust)-int(targetnum))>3:
            print('iter '+str(i))
            if numclust>targetnum:
                resol=resol-0.1
            if numclust<targetnum:  
                resol=resol+0.1
            sc.tl.leiden(adata,resolution=resol,key_added='leiden_1_4')
            numclust=np.max(adata.obs['leiden_1_4'].astype(int))
            i=i+1
    resol=1.4
    sc.tl.louvain(adata,resolution=resol,key_added='louvain_1_4')
    numclust=np.max(adata.obs['louvain_1_4'].astype(int))
    if (numclust-targetnum)>3:
        i=0
        while (int(numclust)-int(targetnum))>3:
            print('iter '+str(i))
            if numclust>targetnum:
                resol=resol-0.1
            if numclust<targetnum:  
                resol=resol+0.1
            sc.tl.leiden(adata,resolution=resol,key_added='louvain_1_4')
            numclust=np.max(adata.obs['louvain_1_4'].astype(int))
            i=i+1
    if default==True:
        return adata,targetnum
    else:
        return adata
    
    
         
def preprocess_adata(adata,save=True,clustering_params={},output_path='output_path'):
     """ Preprocess and cluster the cells in adata given the parameters specified.
   
    Parameters:
    adata (AnnData): AnnData object with the cells of the experiment. 
    save (boolean):whether to save or not the adata object once it has been processed.
    clustering_params(dict): Dictionary where main preprocessing and clustering parameters are inputed
    output_path(str): path where to save the adata object in case that option is selected
     
    

    Returns:
    adata: AnnData object with the preprocessed and clustered cells according to the parameters specified
   """
    
    import matplotlib
    sc.set_figure_params(scanpy=True, dpi=150,figsize=(10,10))
    plt.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    plot_path=output_path+r'/figures/'
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)
    plot_cell_counts(adata,plot_path=plot_path,clustering_params=clustering_params)
    sc.pp.filter_cells(adata,min_counts=clustering_params['min_counts_x_cell'])
    sc.pp.filter_cells(adata,min_genes=clustering_params['min_genes_x_cell'])
    adata.raw=adata
    adata.layers['raw']=adata.X
    sc.pp.normalize_total(adata, target_sum=clustering_params['normalization_target_sum'])
    sc.pp.log1p(adata)
    plt.rcParams['figure.facecolor'] = 'white'
    sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)
    sc.pl.highly_variable_genes(adata,show=False)
    plt.savefig(plot_path+'hvg.pdf')
    if clustering_params['scale']==True:
        sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    plt.rcParams['figure.facecolor'] = 'white'
    sc.pl.pca_variance_ratio(adata, log=True,show=False)
    plt.savefig(plot_path+'pca_variance_ratio.pdf')
    sc.pp.neighbors(adata, n_neighbors=clustering_params['n_neighbors'], n_pcs=clustering_params['n_pcs'])
    if clustering_params['clustering_alg']=='louvain':
        for r in clustering_params['resolutions']:
            sc.tl.louvain(adata,resolution=r,key_added=clustering_params['clustering_alg']+'_'+str(r))
    if clustering_params['clustering_alg']=='leiden':
        for r in clustering_params['resolutions']:
            sc.tl.leiden(adata,resolution=r,key_added=clustering_params['clustering_alg']+'_'+str(r))
    sc.tl.umap(adata,min_dist=clustering_params['umap_min_dist'])
    for r in clustering_params['resolutions']:
        plt.rcParams['figure.facecolor'] = 'white'
        sc.pl.umap(adata,color=[clustering_params['clustering_alg']+'_'+str(r)],size=1,legend_fontsize=8,legend_fontoutline=1,show=False,frameon=False)
        plt.savefig(plot_path+'umap_'+str(r)+'.pdf')   
    sc.pl.umap(adata,color=['sample'],size=1,legend_fontsize=8,legend_fontoutline=1,show=False,frameon=False)
    plt.savefig(plot_path+'umap_sample.pdf')
    for r in clustering_params['resolutions']:
        plt.rcParams['figure.facecolor'] = 'white'
        sc.tl.rank_genes_groups(adata, groupby=clustering_params['clustering_alg']+'_'+str(r), method='wilcoxon',key_added=clustering_params['clustering_alg']+'_'+str(r))
        sc.tl.dendrogram(adata,groupby=clustering_params['clustering_alg']+'_'+str(r))
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, swap_axes=True,show=False,key=clustering_params['clustering_alg']+'_'+str(r))
        plt.savefig(plot_path+'deg_dotplot_'+str(r)+'.pdf')
        sc.pl.rank_genes_groups(adata,n_genes=5,show=False,key=clustering_params['clustering_alg']+'_'+str(r),fontsize=15)
        plt.savefig(plot_path+'deg_'+str(r)+'.pdf')
    for s in adata.obs['sample'].unique():
        adatasub=adata[adata.obs['sample']==s]
        for r in clustering_params['resolutions']:
            sc.pl.spatial(adatasub,color=clustering_params['clustering_alg']+'_'+str(r),spot_size=40)
            plt.savefig(plot_path+'spatial_map_'+str(s)+'_'+str(r)+'.pdf')
    if save==True:
        adata.write(output_path+'combined_processed.h5ad')
    return adata
