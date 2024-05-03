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

def dist_nuc(reads_ctdsub):
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    allds=[]
    for g,n in reads_ctdsub.groupby('cell_id'):
        try:
            hull = ConvexHull(np.array(n.loc[:,['x_location','y_location']]))
            allds.append(np.mean(n.iloc[hull.vertices]['distance']))
        except:
            print()
    return np.median(allds)

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def distance_calc (x1,y1,x2,y2):
    import math
    distance = math.sqrt( ((x1-x2)**2)+((y1-y2)**2) )
    return(distance)

# DISPERSION CALCULATES THE DISTANCE BETWEEN EACH READ AND ITS ASSIGNED CELL
def dispersion(reads_original,adata1):
    reads_assigned=reads_original
    cells_metadata=adata1.obs
    cells_metadata_filt=cells_metadata.loc[cells_metadata['cell_id'].isin(reads_assigned['cell_id']),:]
    def distance_calc (x1,y1,x2,y2):
        import math
        distance = math.sqrt( ((x1-x2)**2)+((y1-y2)**2) )
        return(distance)
    reads_assigned=reads_assigned[reads_assigned['cell_id'].isin(cells_metadata_filt['cell_id'])]


    dictx=dict(zip(cells_metadata_filt['cell_id'],cells_metadata_filt['x_centroid']))
    dicty=dict(zip(cells_metadata_filt['cell_id'],cells_metadata_filt['y_centroid']))
    reads_assigned['x_cell']=list(reads_assigned['cell_id'].map(dictx))
    reads_assigned['y_cell']=list(reads_assigned['cell_id'].map(dicty))
    reads_assigned['distance']=list(np.sqrt(((reads_assigned['x_location']-reads_assigned['x_cell'])**2)+((reads_assigned['y_location']-reads_assigned['y_cell'])**2)))
    return reads_assigned

def entropy(clustering):
    _, counts = np.unique(clustering, return_counts=True)
    proportions = counts / len(clustering)
    return -np.sum(proportions * np.log(proportions))

def compute_vi(ground_truth, predicted):
    mi = mutual_info_score(ground_truth, predicted)
    h_gt = entropy(ground_truth)
    h_pred = entropy(predicted)
    vi_score = h_gt + h_pred - 2 * mi
    return vi_score

def compute_fmi(ground_truth, predicted):
    fmi_score = fowlkes_mallows_score(ground_truth, predicted)
    return fmi_score

from sklearn.metrics import normalized_mutual_info_score
def compute_nmi(ground_truth, predicted):
    nmi_score = normalized_mutual_info_score(ground_truth, predicted)
    return nmi_score

from shapely.geometry import Point, Polygon
def domainassign(plsin,adatadom):
    adatadom.obs['region_annotation']='None'
    plt.figure()
    for sel in plsin['region_annotation'].unique():
        plsub=plsin[plsin['region_annotation']==sel]
        if plsub.shape[0]>2:
    #        plt.figure()
            coord = np.array(plsub[['y','x']]).tolist()
            coord.append(coord[0])
            poli=Polygon(coord)
            xs, ys = zip(*coord) #create lists of x and y values
            plt.plot(xs,ys,color='black') 
            plt.title(sel)
            for n in adatadom.obs.index:
                pnt=Point(adatadom.obs.loc[n,'y_centroid'],adatadom.obs.loc[n,'x_centroid'])
                if pnt.within(poli)==True:
                    adatadom.obs.loc[n,'region_annotation']=sel
    plt.scatter(adatadom.obs['y_centroid'],adatadom.obs['x_centroid'],s=0.5)
    plt.show() # if you need...


def negative_marker_purity_coexpression(adata_sp: AnnData, adata_sc: AnnData, key: str='celltype', pipeline_output: bool=True,minexp:float =0.0):
    """ Negative marker purity aims to measure read leakeage between cells in spatial datasets. 
    
    For this, we calculate the increase in reads assigned in spatial datasets to pairs of genes-celltyes with no/very low expression in scRNAseq
    
    Parameters
    ----------
    adata_sp : AnnData
        Annotated ``AnnData`` object with counts from spatial data
    adata_sc : AnnData
        Annotated ``AnnData`` object with counts scRNAseq data
    key : str
        Celltype key in adata_sp.obs and adata_sc.obs
    pipeline_output : float, optional
        Boolean for whether to use the function in the pipeline or not
    Returns
    -------
    negative marker purity : float
       Increase in proportion of reads assigned in spatial data to pairs of genes-celltyes with no/very low expression in scRNAseq
    """   
    
    # Set threshold parameters
    min_number_cells=10 #minimum number of cells belonging to a cluster to consider it in the analysis
    minimum_exp=0.05 #maximum relative expression allowed in a gene in a cluster to consider the gene-celltype pair the analysis 
    
    # Subset adata_sc to genes of spatial data
    adata_sp = adata_sp[:,adata_sp.var_names.isin(adata_sc.var_names)]
    adata_sc = adata_sc[:,adata_sp.var_names]
    # TMP fix for sparse matrices, ideally we don't convert, and instead have calculations for sparse/non-sparse
    #for a in [adata_sc, adata_sp]:
    #    if issparse(a.X):
    #        a.X = a.X.toarray()
    print(adata_sp.shape)
    print(adata_sc.shape)
    
    print(adata_sc.var.index)
    # Get mean expression per cell type
    try:
        exp_sc=pd.DataFrame(adata_sc.X.todense(),columns=adata_sc.var_names)
    except:
        exp_sc=pd.DataFrame(adata_sc.X,columns=adata_sc.var_names)
    try: 
        exp_sp = pd.DataFrame(adata_sp.X.todense(),columns=adata_sp.var_names)
    except:
        exp_sp = pd.DataFrame(adata_sp.X,columns=adata_sp.var_names)
    print(exp_sc.shape)
    #exp_sc['celltype'] = list(adata_sc.obs[key])
    #exp_sp['celltype'] = list(adata_sp.obs[key])
    # instead of measuring the purity by cell type, we check it based on the proportion of posive cells
    
    #mean_celltype_sc = exp_sc.groupby('celltype').mean()
    #mean_celltype_sp = exp_sp.groupby('celltype').mean()
    mean_celltype_sc=coexpression_calculation(exp_sc,min_exp=minexp)
    mean_celltype_sp=coexpression_calculation(exp_sp,min_exp=minexp)
    
    # Get mean expressions relative to mean over cell types (this will be used for filtering based on minimum_exp)
    mean_ct_sc_rel = mean_celltype_sc#.div(mean_celltype_sc.mean(axis=0),axis=1)
    mean_ct_sp_rel = mean_celltype_sp#.div(mean_celltype_sp.mean(axis=0),axis=1)
    
    # Get normalized mean expressions over cell types (this will be summed up over negative cell types)
    mean_ct_sc_norm = mean_celltype_sc#.div(mean_celltype_sc.sum(axis=0),axis=1)
    mean_ct_sp_norm = mean_celltype_sp#.div(mean_celltype_sp.sum(axis=0),axis=1)
    # Explanation: The idea is to measure which ratio of mean expressions is shifted towards negative clusters.
    #              With this approach the metric stays between 0 and 1
    
    commongenes=mean_ct_sc_rel.index
    # Get gene-cell type pairs with negative marker expression
    neg_marker_mask = np.array(mean_ct_sc_rel < minimum_exp)
    sns.clustermap(mean_ct_sc_rel.astype(float),vmax=0.1)
    # Return nan if no negative markers were found
    if np.sum(neg_marker_mask) < 1:
        print("No negative markers were found in the sc data reference.")
        negative_marker_purity = 'nan'
        if pipeline_output==True:
            return negative_marker_purity
        else:
            return negative_marker_purity, None, None
    
    # Get marker expressions in negative marker-cell type pairs
    lowvals_sc = np.full_like(neg_marker_mask, np.nan, dtype=np.float32)
    lowvals_sc = mean_ct_sc_norm.values[neg_marker_mask]
    lowvals_sp = mean_ct_sp_norm.values[neg_marker_mask]
    
    # Take the mean over the normalized expressions of the genes' negative cell types
    mean_sc_lowexp = np.nanmean(lowvals_sc)
    mean_sp_lowexp = np.nanmean(lowvals_sp)
    
    # Calculate summary metric
    #negative_marker_purity = 1
    #if mean_sp_lowexp > mean_sc_lowexp:
    #    negative_marker_purity -= (mean_sp_lowexp - mean_sc_lowexp)
    lowvals_diff=(lowvals_sp-lowvals_sc)
    lowvals_diff[lowvals_diff<0]=0
    negative_marker_purity=1-np.mean(lowvals_diff)
    
    if pipeline_output:
        return negative_marker_purity
    else:
        # Calculate per gene and per cell type purities
        purities = (mean_ct_sp_norm - mean_ct_sc_norm)#.cip(0,None)
        purities[~neg_marker_mask] = np.nan
        purities = purities.loc[~(purities.isnull().all(axis=1)), ~(purities.isnull().all(axis=0))]
        purity_per_gene = purities.mean(axis=0, skipna=True)
        purity_per_celltype = purities.mean(axis=1, skipna=True)
        return negative_marker_purity, purity_per_gene, purity_per_celltype,lowvals_sc,lowvals_sp,commongenes

    
def coexpression_calculation(exp,min_exp=0):
    coexpression=pd.DataFrame(index=exp.columns,columns=exp.columns)
    for col in tqdm(exp.columns):
        sel=exp.loc[:,col]>min_exp
        positive_cells=exp.loc[sel,:]
        coexpression.loc[:,col]=np.sum(positive_cells>min_exp)/positive_cells.shape[0]
    coexpression=coexpression.fillna(1)
    return coexpression


def alphashape_fun(points,alpha=0.1):
    alpha_shape = alphashape.alphashape(points, alpha)
    area = alpha_shape.area
    try:
        ax.add_patch(plt.PolygonPatch(alpha_shape, alpha=0.2, color='orange', label='Alpha Shape'))
    except:
        ss='ss'

    # Print the area
    print("Area of the Alpha Shape (Concave Hull):", area)
    return area

def svf_moranI(adata1,sample_key='sample',radius=50.0):
    anndata_list=[]
    for sample in adata.obs[sample_key].unique():
        adata_copy_int = adata[adata.obs[sample_key]==sample]
        sq.gr.spatial_neighbors(adata_copy_int,radius=radius,coord_type ='generic')
        anndata_list.append(adata_copy_int) 
    adata1=sc.concat(anndata_list,join='outer',pairwise=True) 
    sq.gr.spatial_autocorr(adata1, mode="moran")
    hs_results=adata1.uns["moranI"]
    hs_results['rank']=list(hs_results['I'].rank())
    hs_results=hs_results.loc[:,['pval_norm','pval_norm_fdr_bh','rank']]
    hs_results.columns=['Pval','FDR','rank']
    return adata1,hs_results