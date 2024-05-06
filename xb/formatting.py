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


def format_xenium_adata(path,tag,output_path):
    """ Format xenium data (output from the machine) to adata format, using the original Xenium format (pre-release)
   
    Parameters:
    path(str): path to the folder where the output of the Xenium machine is stored
    tag(str): sample tag to be added to be added to all cells formated from the section. 
    output_path(str): path where to store the resulting adata object

    Returns:
    adata: AnnData object with the formated cells

   """
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
    adata.obsm['spatial']=np.array(adata.obs.loc[:,['x_centroid','y_centroid']])
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
    if os.path.isfile(path+'/background.tiff')==False:
        format_background(path)
    IM=tf.TiffFile(path+'/background.tiff')
    position1_series = IM.series[0]
    position1_series.axes
    position1 = position1_series.asarray()
    image_downsize_fact=1/(2000/np.max(position1.shape))
    pos1_resized=np.resize(position1,(position1.shape/image_downsize_fact).astype(int))
    adata.uns={"spatial":{tag:{"scalefactors":{"tissue_um_to_pixel":1/0.2125,"tissue_hires_scalef":1/(0.2125*image_downsize_fact)},"images":{"hires":pos1_resized}}}} 
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
        kmeans3=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_3_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans3_clusters']=list(kmeans3['Cluster'].astype(str))
        kmeans4=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_4_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans4_clusters']=list(kmeans4['Cluster'].astype(str))
        kmeans5=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_5_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans5_clusters']=list(kmeans5['Cluster'].astype(str))
        kmeans6=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_6_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans6_clusters']=list(kmeans6['Cluster'].astype(str))
        kmeans7=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_7_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans7_clusters']=list(kmeans7['Cluster'].astype(str))
        kmeans8=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_8_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans8_clusters']=list(kmeans8['Cluster'].astype(str))
        kmeans9=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_9_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans9_clusters']=list(kmeans9['Cluster'].astype(str))
        kmeans10=pd.read_csv(path+'/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv',index_col=0)
        adata.obs['kmeans10_clusters']=list(kmeans10['Cluster'].astype(str))
    except:
        print('UMAP and clusters_could not be recovered')
    adata.obs['cell_id']=adata.obs['cell_id'].astype(str)
    adata.write(output_path+tag+'.h5ad')
    return adata
    

def format_xenium_adata_2023(path,tag,output_path):
    """ Format xenium data (output from the machine) to adata format, considerin the format used by Xenium in Q1 2023
   
    Parameters:
    path(str): path to the folder where the output of the Xenium machine is stored
    tag(str): sample tag to be added to be added to all cells formated from the section. 
    output_path(str): path where to store the resulting adata object

    Returns:
    adata: AnnData object with the formated cells

   """
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
    # Opening JSON file
    f = open(path+'/gene_panel.json')
    # returns JSON object as 
    # a dictionary
    data = json.load(f)
    # Iterating through the json
    # list
    geness=[]
    idss=[]
    descriptorss=[]
    for r in range(len(data['payload']['targets'])):
        geness.append(data['payload']['targets'][r]['type']['data']['name'])
        try:
            idss.append(data['payload']['targets'][r]['type']['data']['id'])
        except:
            idss.append('newid_'+str(r)) 
        try:
            descriptorss.append(data['payload']['targets'][r]['type']['descriptor'])
        except:
            descriptorss.append('other')
    # Closing file
    f.close()

    dict_inpanel=dict(zip(geness,descriptorss))
    dict_ENSEMBL=dict(zip(geness,idss))
    adata.var['Ensembl ID']=adata.var['gene_id'].map(dict_ENSEMBL)
    adata.var['in_panel']=adata.var['gene_id'].map(dict_inpanel)
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
    adata.X=sp.sparse.csr_matrix(adata.X)
    adata.write(output_path+tag+'.h5ad')
    return adata


def format_xenium_adata_mid_2023(path,tag,output_path):
    """ Format xenium data (output from the machine) to adata format, considerin the format used by Xenium at Q2 2023
   
    Parameters:
    path(str): path to the folder where the output of the Xenium machine is stored
    tag(str): sample tag to be added to be added to all cells formated from the section. 
    output_path(str): path where to store the resulting adata object

    Returns:
    adata: AnnData object with the formated cells

   """
    if os.path.isdir(path+'/cell_feature_matrix')==False:
    # Path of the file
        extr=path+'/cell_feature_matrix.tar.gz'
        # Target directory
        extract_dir = path
        # Unzip the file 
        shutil.unpack_archive(extr, extract_dir)
        print('First decompressing done')
    
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
    f = open(path+'/gene_panel.json') 
    data = json.load(f)
    geness=[]
    idss=[]
    descriptorss=[]
    for r in range(len(data['payload']['targets'])):
        geness.append(data['payload']['targets'][r]['type']['data']['name'])
        try:
            idss.append(data['payload']['targets'][r]['type']['data']['id'])
        except:
            idss.append('newid_'+str(r)) 
        try:
            descriptorss.append(data['payload']['targets'][r]['type']['descriptor'])
        except:
            descriptorss.append('other')
    # Closing file
    f.close()

    dict_inpanel=dict(zip(geness,descriptorss))
    dict_ENSEMBL=dict(zip(geness,idss))
    adata.var['Ensembl ID']=adata.var['gene_id'].map(dict_ENSEMBL)
    adata.var['in_panel']=adata.var['gene_id'].map(dict_inpanel)
    transcripts=pd.read_csv(path+'/transcripts.csv',index_col=0)
    adata.uns['spots']=transcripts
    if os.path.isdir(path+'/analysis')==False:
    # Path of the file
        extr=path+'/analysis.tar.gz'
        # Target directory
        extract_dir = path
        # Unzip the file 
        shutil.unpack_archive(extr, extract_dir)
        print('Analysis files decompressed')
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
    adata.X=sp.sparse.csr_matrix(adata.X)
    adata.uns['spots']=adata.uns['spots'].fillna(0)
    adata.uns['spots']['fov_name']=adata.uns['spots']['fov_name'].astype(str)
    adata.write(output_path+tag+'.h5ad')
    return adata


def format_background(path):
    """ Format OME-TIFF background mipped image to .tiff image
   
    Parameters:
    path(str): path to the folder where the output of the Xenium machine is stored
   
    Returns:
    None

   """
    IM=tf.TiffFile(path+'/morphology_mip.ome.tif')
    position1_series = IM.series[0]
    position1_series.axes
    position1 = position1_series.asarray()
    tf.imwrite(path+'/background.tiff',position1)
    
    
def keep_nuclei(adata1,overlaps_nucleus=1):
    """ Redefine cells in AnnData to keep only nuclear reads 
   
    Parameters:
    adata1(AnnData): AnnData object with the cells of the experiment.
    overlaps_nucleus(int): whether to keep only nuclear reads only (1) or cytoplasmic reads (0) in the redefinition of cells 
    
    Returns:
    adata: AnnData object with the formated cells

   """
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


def generate_random_color_variation(base_color, deviation=0.17):
     """ Generate variations of a reference color
   
    Parameters:
    base_color (str):reference hex color
    deviation(float): deviation from the base color that the resulting color should have.
   
    Returns:
    modified_hex_color(str):resulting hex color

   """
    
    base_rgb = mcolors.hex2color(base_color)
    h, s, v = mcolors.rgb_to_hsv(base_rgb)

    # Make random adjustments to the hue, saturation, and value
    h += random.uniform(-deviation, deviation)
    s = max(0, min(1, s + random.uniform(-deviation, deviation)))
    v = max(0, min(1, v + random.uniform(-deviation, deviation)))

    # Ensure the values are within the valid HSV range
    h = h % 1.0

    modified_rgb = mcolors.hsv_to_rgb((h, s, v))
    modified_hex_color = mcolors.to_hex(modified_rgb)

    return modified_hex_color


def format_data_neighs(adata,sname,condit,neighs=10):
    """ Redefine the expression of cells in adata by counting the neighnoring cell types of each cell

    Parameters:
    adata (AnnData): AnnData object with the cells of the experiment
    sname(str): column in adata.obs where the cluster assigned to each cells are stored
    neighs(int): number of neighbors to consider when computing neighboring cells
    
    Returns:
    adata1 (AnnData): AnnData object with neighboring cell types included in a cell-by-celltype matrix
    
   """
    
    try:
        adata.obsm['spatial']
    except:
        adata.obsm["spatial"]=np.array([adata.obs.X,adata.obs.Y]).transpose().astype('float64')
    adata_copy_int=adata
    sq.gr.spatial_neighbors(adata_copy_int,n_neighs=neighs)
    result=np.zeros([adata.shape[0],len(adata_copy_int.obs[sname].unique())])
    n=0
    tr=adata_copy_int.obsp['spatial_distances'].transpose()
    tr2=tr>0
    from tqdm import tqdm
    for g in tqdm(adata_copy_int.obs[sname].unique()):
        epv=adata_copy_int.obs[sname]==g*1
        opv=list(epv*1)
        result[:,n]=tr2.dot(opv)
        n=n+1
    expmat=pd.DataFrame(result,columns=adata_copy_int.obs[sname].unique())
    adata1=sc.AnnData(expmat,obs=adata.obs)
    #adata1.obs['sample']=condit
    adata1.obs['condition']=condit
    return adata1

def format_data_neighs_colapse(adata,sname,condit,neighs=10):
     """ Redefine the expression of cells in adata by collapsing the expression of its neighbors into each cell (a.k.a pseudobining)

    Parameters:
    adata (AnnData): AnnData object with the cells of the experiment
    sname(str): column in adata.obs where sample is stored
    condit(str): column in adata.obs where the sample each cell belongs to is stored
    neighs(int): number of neighbors to consider when collapsing the expression of neighboring cells
    
    Returns:
    adata1 (AnnData): AnnData object with expression of cells collapsed from neighboring cells
    
   """
    adata.obsm["spatial"]=np.array([adata.obs.X,adata.obs.Y]).transpose().astype('float64')
    adata_copy_int=adata
    sq.gr.spatial_neighbors(adata_copy_int,n_neighs=neighs)
    result=np.zeros([adata.shape[0],adata.shape[1]])
    n=0
    tr=adata_copy_int.obsp['spatial_distances'].transpose()
    tr2=tr>0
    exp=adata_copy_int.to_df()
    from tqdm import tqdm
    #tdd=tr2.todense()
    for i in tqdm(range(0,adata_copy_int.to_df().shape[0])):
        result[i,:]=np.sum(exp[tr2[i,:].todense().transpose()],axis=0)
    adata1=sc.AnnData(result,obs=adata.obs,var=adata.var)
    return adata1



def format_xenium_adata_final(path,tag,output_path,use_parquet=True,save=True):
    """ Format xenium data (output from the machine) to adata format using the official up-to-date Xenium format.
   
    Parameters:
    path(str): path to the folder where the output of the Xenium machine is stored, if requested
    tag(str): sample tag to be added to be added to all cells formated from the section. 
    output_path(str): path where to store the resulting adata object
    use_parquet(boolean): whether to use parquet files as an input to generate the AnnData File. (it's way faster)
    save(boolean): whether to save the resulting object

    Returns:
    adata: AnnData object with the formated cells

   """
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
    f = open(path+'/gene_panel.json')
    data = json.load(f)
    geness=[]
    idss=[]
    descriptorss=[]
    for r in range(len(data['payload']['targets'])):
        geness.append(data['payload']['targets'][r]['type']['data']['name'])
        try:
            idss.append(data['payload']['targets'][r]['type']['data']['id'])
        except:
            idss.append('newid_'+str(r)) 
        try:
            descriptorss.append(data['payload']['targets'][r]['type']['descriptor'])
        except:
            descriptorss.append('other')
    f.close()

    dict_inpanel=dict(zip(geness,descriptorss))
    dict_ENSEMBL=dict(zip(geness,idss))
    adata.var['Ensembl ID']=adata.var['gene_id'].map(dict_ENSEMBL)
    adata.var['in_panel']=adata.var['gene_id'].map(dict_inpanel)
    if use_parquet==True:
        transcripts=pd.read_parquet(path+'/transcripts.parquet')
    else:
        if os.path.isfile(path+'/transcripts.csv')==False:
            with gzip.open(path+'/transcripts.csv.gz', 'rb') as f_in:
                with open(path+'/transcripts.csv', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        transcripts=pd.read_csv(path+'/transcripts.csv')
    adata.uns['spots']=transcripts
    
    UMAP=pd.read_csv(path+'/analysis/umap/gene_expression_2_components/projection.csv',index_col=0)
    id2UMAP1=dict(zip(UMAP.index,UMAP.iloc[:,0]))
    id2UMAP2=dict(zip(UMAP.index,UMAP.iloc[:,1]))
    adata.obsm['X_umap_original']=np.array([adata.obs['cell_id'].map(id2UMAP1),adata.obs['cell_id'].map(id2UMAP2)]).transpose()
    PCA=pd.read_csv(path+'/analysis/pca/gene_expression_10_components/projection.csv',index_col=0)
    id2PCA1=dict(zip(PCA.index,PCA.iloc[:,0]))
    id2PCA2=dict(zip(PCA.index,PCA.iloc[:,1]))
    adata.obsm['X_pca_10_comp']=np.array([adata.obs['cell_id'].map(id2PCA1),adata.obs['cell_id'].map(id2PCA2)]).transpose()
    dici={'/analysis/clustering/gene_expression_graphclust/clusters.csv':'graph_clusters',
           '/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv':'kmeans2',
           '/analysis/clustering/gene_expression_kmeans_3_clusters/clusters.csv':'kmeans3',
           '/analysis/clustering/gene_expression_kmeans_4_clusters/clusters.csv':'kmeans4',
           '/analysis/clustering/gene_expression_kmeans_5_clusters/clusters.csv':'kmeans5',
           '/analysis/clustering/gene_expression_kmeans_6_clusters/clusters.csv':'kmeans6',
           '/analysis/clustering/gene_expression_kmeans_7_clusters/clusters.csv':'kmeans7',
           '/analysis/clustering/gene_expression_kmeans_8_clusters/clusters.csv':'kmeans8',
           '/analysis/clustering/gene_expression_kmeans_9_clusters/clusters.csv':'kmeans9',
          '/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv':'kmeans10'
    }
    for subpath in dici.keys():
        clust=pd.read_csv(path+subpath,index_col=0)
        clustdict=dict(zip(clust.index,clust.iloc[:,0]))
        adata.obs[dici[subpath]]=adata.obs['cell_id'].map(clustdict).astype(str)
    adata.X=sp.sparse.csr_matrix(adata.X)
    adata.uns['spots']=adata.uns['spots'].fillna(0)
    adata.uns['spots']['fov_name']=adata.uns['spots']['fov_name'].astype(str)
    adata.obsm['spatial']=np.array([adata.obs['x_centroid'],adata.obs['y_centroid']]).transpose()
    adata.obs['sample']=tag
    if save==True:
        adata.write(output_path+tag+'_original.h5ad')
    return adata

def keep_nuclei_and_quality(adata1,tag:str,max_nucleus_distance=1,min_quality=20,save=True,output_path=''):
     """ Redefine cell expression based on nuclei expression an quality of detected reads
   
    Parameters:
    adata1 (AnnData): AnnData object with the cells of the experiment before filtereing reads based on quality or nuclear/non-nuclear
    tag (str): sample tag to added in the name of the saved filed, if needed
    save(boolean): whether to save the resulting files
    output_path(str): if needed, where to save the resulting files
    max_nucleus_distance(float): Maximum distance from the nuclei for reads to be kept in redefined cells
    min_quality(float): Define minimum quality (qv) of reads to keep in the analysis 
   
    Returns:
    adata1nuc(AnnData): AnnData object with the cells redefined based to input parameters

   """
    if max_nucleus_distance==0:
        subset1=adata1.uns['spots'].loc[adata1.uns['spots']['overlaps_nucleus']==overlaps_nucleus,:]
    if max_nucleus_distance>0:
        subset1=adata1.uns['spots'].loc[adata1.uns['spots']['nucleus_distance']<max_nucleus_distance,:]
    subset1=subset1[subset1['qv']>min_quality]
    ct1=pd.crosstab(subset1['cell_id'],subset1['feature_name'])
    adataobs=adata1.obs.loc[adata1.obs['cell_id'].isin(ct1.index),:]
    av=adata1.var[adata1.var['gene_id'].isin(ct1.columns)]#.isin(adata1.var.index)
    ct1=ct1.loc[:,av['gene_id']]
    adataobs.index=adataobs['cell_id']
    adataobs.index.name='ind'
    ct1=ct1.loc[ct1.index.isin(adataobs['cell_id']),:]
    adata1nuc=sc.AnnData(np.array(ct1),obs=adataobs,var=av)
    if save==True:
        adata1nuc.write(output_path+tag+'_filtered.h5ad')
    adata1nuc.obsm['spatial']=np.array([adata1nuc.obs['x_centroid'],adata1nuc.obs['y_centroid']]).transpose()
    return adata1nuc

def format_to_adata(files:list,output_path:str,use_parquet=True,save=False,max_nucleus_distance=0,min_quality=10):
    """ Format xenium datasets (outputs from the machine, up to date 2024) to adata files and filter reads based on quality parameters
   
    Parameters:
    files(list): list including the paths where  the Xenium outputs are saved for each sample (output from the machine)
    output_path(str): path where to store the resulting adata object
    use_parquet(boolean): whether to use parquet files as an input to generate the AnnData File. (it's way faster)
    save(boolean): whether to save the resulting object
    max_nucleus_distance: Maximum distance from the nuclei for reads to be kept in redefined cells
    min_quality(float): Define minimum quality (qv) of reads to keep in the analysis 

    Returns:
    adata: AnnData object with the formated cells with only reads that passed the filters established

   """
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    alladata=[]
    for f in files:
        tag=f.split('/')[-1]
        print(f'Formatting {tag}')
        adata=format_xenium_adata_final(f,tag,output_path,use_parquet=use_parquet,save=save)
        print('Filter reads')
        adata_f1=keep_nuclei_and_quality(adata,tag=tag,max_nucleus_distance=max_nucleus_distance,min_quality=min_quality,save=save,output_path=output_path)
        #xf.format_background(f)
        alladata.append(adata_f1)
    adata=sc.concat(alladata,join='outer')
    adata.var=alladata[0].var
    adata.obs['total_counts']=np.sum(adata.X,axis=1)
    adata.obs['expressed_genes']=np.sum(adata.X>0,axis=1)
    adata.obs['unique_cell_id']=adata.obs['cell_id'].astype(str)+'_'+adata.obs['sample'].astype(str)
    adata.var['ENSMBL_ID']=adata.var.index
    adata.var.index=adata.var['gene_id']
    if save==True:
        adata.write(output_path+'combined_filtered.h5ad')
    return adata
