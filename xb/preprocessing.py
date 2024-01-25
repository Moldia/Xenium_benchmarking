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

def main_preprocessing(adata,target_sum=100,mincounts=10,mingenes=3,neigh=15,npc=0,nuc=1,scale=False,hvg=False,default=False,total_clusters=30,norm=True,lg=True):
#    print(adata.shape)
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
#    sc.tl.leiden(adata,resolution=2.2,key_added='leiden_2_2')
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
    