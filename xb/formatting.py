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


def format_xenium_adata(path,tag,output_path):
    #decompress
    if os.path.isfile(path+'/transcripts.csv')==False:
        with gzip.open(path+'/transcripts.csv.gz', 'rb') as f_in:
            with open(path+'/transcripts.csv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    if os.path.isfile(path+'/cell_feature_matrix/barcodes.tsv')==False:
        with gzip.open(path+'/cell_feature_matrix/barcodes.tsv.gz', 'rb') as f_in:
            with open(path+'/cell_feature_matrix/barcodes.tsv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    if os.path.isfile(path+'/cell_feature_matrix/features.tsv')==False:
        with gzip.open(path+'/cell_feature_matrix/features.tsv.gz', 'rb') as f_in:
            with open(path+'/cell_feature_matrix/features.tsv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    if os.path.isfile(path+'/cell_feature_matrix/matrix.mtx')==False:
        with gzip.open(path+'/cell_feature_matrix/matrix.mtx.gz', 'rb') as f_in:
            with open(path+'/cell_feature_matrix/matrix.mtx', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    if os.path.isfile(path+'/cells.csv')==False:
        with gzip.open(path+'/cells.csv.gz', 'rb') as f_in:
            with open(path+'/cells.csv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    a = mmread(path+'/cell_feature_matrix/matrix.mtx')
    ad=a.todense()
    cell_info=pd.read_csv(path+r"/cells.csv")
    features=pd.read_csv(path+'/cell_feature_matrix/features.tsv',sep='\t',header=None,index_col=0)
    barcodes=pd.read_csv(path+'/cell_feature_matrix/barcodes.tsv',header=None,index_col=0)
    adata=sc.AnnData(ad.transpose(),obs=cell_info,var=features)
    adata.var.index.name='index'
    adata.var.columns=['gene_id','reason_of_inclusion']
    panel_info=pd.read_csv(path+'/panel.tsv',sep='\t')
    try:
        panel_info['Gene']
    except:
        panel_info['Gene']=panel_info['Name']
    dict_annotation=dict(zip(panel_info['Gene'],panel_info['Annotation']))
    dict_ENSEMBL=dict(zip(panel_info['Gene'],panel_info['Ensembl ID']))
    adata.var['Annotation']=adata.var.index.map(dict_annotation)
    adata.var['Ensembl ID']=adata.var.index.map(dict_ENSEMBL)
    adata.var['in_panel']=adata.var.index.isin(panel_info['Gene'])
    transcripts=pd.read_csv(path+'/transcripts.csv',index_col=0)
    adata.uns['spots']=transcripts
    try:
        UMAP=pd.read_csv(path+'/analysis/umap/gene_expression_2_components/projection.csv',index_col=0)
        adata.obsm['X_umap']=np.array(UMAP)
        TSNE=pd.read_csv(path+'/analysis/tsne/gene_expression_2_components/projection.csv',index_col=0)
        adata.obsm['X_tsne']=np.array(TSNE)
        PCA=pd.read_csv(path+'/analysis/PCA/gene_expression_10_components/projection.csv',index_col=0)
        adata.obsm['X_pca']=np.array(PCA)
        clusters=pd.read_csv(path+'/analysis/clustering/gene_expression_graphclust/clusters.csv',index_col=0)
        adata.obs['graph_clusters']=list(clusters['Cluster'].astype(str))
        kmeans2=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans2_clusters']=list(kmeans2['Cluster'].astype(str))
        kmeans3=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans3_clusters']=list(kmeans3['Cluster'].astype(str))
        kmeans4=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans4_clusters']=list(kmeans4['Cluster'].astype(str))
        kmeans5=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans5_clusters']=list(kmeans5['Cluster'].astype(str))
        kmeans6=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans6_clusters']=list(kmeans6['Cluster'].astype(str))
        kmeans7=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans7_clusters']=list(kmeans7['Cluster'].astype(str))
        kmeans8=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans8_clusters']=list(kmeans8['Cluster'].astype(str))
        kmeans9=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans9_clusters']=list(kmeans9['Cluster'].astype(str))
        kmeans10=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans10_clusters']=list(kmeans10['Cluster'].astype(str))
    except:
        print('UMAP and clusters_could not be recovered')
    adata.write(output_path+tag+'.h5ad')
    return adata
    

def format_background(path):
    IM=tf.TiffFile(path+'/morphology_mip.ome.tif')
    position1_series = IM.series[0]
    position1_series.axes
    position1 = position1_series.asarray()
    tf.imwrite(path+'/background.tiff',position1)
    
    
def keep_nuclei(adata1,overlaps_nucleus=1):
    subset1=adata1.uns['spots'].loc[adata1.uns['spots']['overlaps_nucleus']==overlaps_nucleus,:]
    ct1=pd.crosstab(subset1['cell_id'],subset1['feature_name'])
    adataobs=adata1.obs.loc[adata1.obs['cell_id'].isin(ct1.index),:]
    ct1=ct1.loc[:,adata1.var.index]
    adataobs.index=adataobs['cell_id']
    adataobs.index.name='ind'
    ct1=ct1.loc[ct1.index.isin(adataobs['cell_id']),:]
    adata1nuc=sc.AnnData(np.array(ct1),obs=adataobs,var=adata1.var)
    return adata1nuc
    
def cell_area(adata_sp: AnnData,pipeline_output=True):
    """Calculates the area of the region imaged using convex hull and divide total number of cells/area. XY position should be in um2"
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    density : float
       Cell density (cells/um)
    """   
    hull = ConvexHull(np.array(adata_sp.uns['spots'].loc[:,['x','y']]))
    area=hull.area#*10e6
    return area