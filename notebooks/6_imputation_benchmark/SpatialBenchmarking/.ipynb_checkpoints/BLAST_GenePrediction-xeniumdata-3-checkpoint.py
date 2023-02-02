#!/usr/bin/env python
# coding: utf-8

# ### These jupyter can reproduced the results of figure2 in our papers.
# ### you can use these functions we defined to evaluate the performance of each method with 10-X cross-validation.
# ### If you want to use each method to analyse youe own data, please see the Tutorial.ipynb

# # 

# ### Please note: We used raw count as input for predictions in this ipynb !!!
# 
# ### In our manuscripts, we tested four schemes of input expression matrices: 
# ### (1) raw expression matrix of spatial data and raw expression matrix of scRNA-seq data (R-R); 
# ### (2) normalized expression matrix of spatial data and raw expression matrix of scRNA-seq data (N-R); 
# ### (3) raw expression matrix of spatial data and normalized expression matrix of scRNA-seq data (R-N); 
# ### (4) normalized expression matrix of spatial data and normalized expression matrix of scRNA-seq data (N-N).
# 
# ### For N-N and N-R input, we calculate the accuracy by comapring the prediction result with normalized matrix of spatial data.
# ### For R-R and R-N input, we calculate the accuracy by comapring the prediction result with raw expression matrix of spatial data.

# In[22]:


import numpy as np
import pandas as pd
import sys
import pickle
import os
import time as tm
from functools import partial
import scipy.stats as st
from scipy.stats import wasserstein_distance
import scipy.stats
import copy
from sklearn.model_selection import KFold
import pandas as pd
import multiprocessing
import matplotlib as mpl 
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

from scipy.spatial import distance_matrix
from sklearn.metrics import matthews_corrcoef
from scipy import stats
import seaborn as sns

from scipy.spatial.distance import cdist
import h5py
from scipy.stats import spearmanr

import time
import sys
import tangram as tg
from os.path import join
from IPython.display import display


# In[23]:


def SpaGE_impute(K):
    print ('We run SpaGE for this data\n')
    sys.path.append("Extenrnal/SpaGE-master/")
    from SpaGE.main import SpaGE
    global RNA_data, Spatial_data, train_gene, predict_gene, iterations
    
    RNA_data = RNA_data.loc[(RNA_data.sum(axis=1) != 0)] # removes genes that are in the train list, causing error
    RNA_data = RNA_data.loc[(RNA_data.var(axis=1) != 0)]
        
#     train_gene = []
#     predict_gene = []

#     for i in range(iterations):
#         train_gene.append(list(RNA_data.index&Spatial_data.columns))
#         predict_gene.append(list(set(RNA_data.index) - set(Spatial_data.columns))[:20])

    # train = train_gene[K]
    # predict = predict_gene[K]
    
    train = np.array(train_gene[K])
    predict = np.array(predict_gene[K])
        
    print('train gene')

    pv = len(train)/2
    Spatial = Spatial_data[train]
    Img_Genes = SpaGE(Spatial, RNA_data.T, n_pv = int(pv), genes_to_predict = predict)
    
    print('Img_Genes')
    
    result = Img_Genes[predict]
    return result

def gimVI_impute(K): # from tutorial
    print ('We run gimVI for this data\n')
    import scvi
    import scanpy as sc
    from scvi.model import GIMVI
    import torch
    from torch.nn.functional import softmax, cosine_similarity, sigmoid
    global RNA_data_adata, Spatial_data_adata, train_gene, predict_gene
    
    test_list = predict_gene[K]
    train_list = train_gene[K]
    
    print('train list')

    Genes  = train_list.copy()
    Genes.extend(test_list)
    rand_test_gene_idx = [Genes.index(x) for x in test_list]
    n_genes = len(Genes)
    rand_train_gene_idx = [Genes.index(x) for x in train_list]
    rand_train_genes = np.array(Genes)[rand_train_gene_idx]
    rand_test_genes = np.array(Genes)[rand_test_gene_idx]
    
    spatial_data_partial = Spatial_data_adata[:, rand_train_genes]
        
    sc.pp.filter_cells(spatial_data_partial, min_counts= 0)

    seq_data = copy.deepcopy(RNA_data_adata)

    seq_data = seq_data[:, Genes]
    sc.pp.filter_cells(seq_data, min_counts = 0)
    scvi.data.setup_anndata(spatial_data_partial)
    scvi.data.setup_anndata(seq_data)

    model = GIMVI(seq_data, spatial_data_partial)
    model.train(200)

    _, imputation = model.get_imputed_values(normalized=False)
    imputed = imputation[:, rand_test_gene_idx]
    result = pd.DataFrame(imputed, columns=rand_test_genes)

    print (test_list ==  rand_test_genes)
    print (train_list ==  rand_train_genes)
    
    return result
                       
def novoSpaRc_impute(K): # modified from tutorial
    print ('We run novoSpaRc for this data\n')
    import novosparc as nc
    global RNA_data, Spatial_data, locations, train_gene, predict_gene
    
    test_list = predict_gene[K]
    train_list = train_gene[K]
    
    locations = np.array(locations[['x_centroid','y_centroid']])
    
    gene_names = np.array(RNA_data.index.values)
    dge = RNA_data.values
    dge = dge.T
    num_cells = dge.shape[0]
    print ('number of cells and genes in the matrix:', dge.shape)

    hvg = np.argsort(np.divide(np.var(dge,axis=0),np.mean(dge,axis=0)+0.0001))
    dge_hvg = dge[:,hvg[-2000:]]
        
    num_locations = locations.shape[0]
    
    print('num locations: ', num_locations)
        
    p_location, p_expression = nc.rc.create_space_distributions(num_locations, num_cells)
    cost_expression, cost_locations = nc.rc.setup_for_OT_reconstruction(dge_hvg,locations,num_neighbors_source = 5,num_neighbors_target = 5)

    insitu_matrix = np.array(Spatial_data[train_list])
    insitu_genes = np.array(Spatial_data[train_list].columns)
    test_genes = np.array(test_list)
        
    markers_in_sc = np.array([], dtype='int')
    for marker in insitu_genes:
        marker_index = np.where(gene_names == marker)[0]
        if len(marker_index) > 0:
            markers_in_sc = np.append(markers_in_sc, marker_index[0])

    cost_marker_genes = cdist(dge[:, markers_in_sc]/np.amax(dge[:, markers_in_sc]),insitu_matrix/np.amax(insitu_matrix))
    alpha_linear = 0.5
        
    gw = nc.rc._GWadjusted.gromov_wasserstein_adjusted_norm(cost_marker_genes, cost_expression, cost_locations,alpha_linear, p_expression, p_location,'square_loss', epsilon=5e-3, verbose=True)

    sdge = np.dot(dge.T, gw)
    imputed = pd.DataFrame(sdge,index=RNA_data.index)
    result = imputed.loc[test_genes]
    result = result.T
    
    return result
   
                       
def SpaOTsc_impute(K):
    print ('We run SpaOTsc for this data\n')
    from spaotsc import SpaOTsc
    global RNA_data, Spatial_data, locations, train_gene, predict_gene
    
    test_list = predict_gene[K]
    train_list = train_gene[K]

    print('train list')
    
    locations = np.array(locations[['x_centroid','y_centroid']])
                          
    df_sc = RNA_data.T
    df_IS = Spatial_data
    pts = locations
    is_dmat = distance_matrix(pts, pts)
    df_is = df_IS.loc[:,train_list]
    
    print('data')
    
    gene_is = df_is.columns.tolist()
    gene_sc = df_sc.columns.tolist()
    gene_overloap = list(set(gene_is).intersection(gene_sc))
    a = df_is[gene_overloap]
    b = df_sc[gene_overloap]
    
    print('gene')
    
    rho, pval = stats.spearmanr(a, b,axis=1)
    rho[np.isnan(rho)]=0
    mcc=rho[-(len(df_sc)):,0:len(df_is)]
    C = np.exp(1 - mcc)
        
    issc = SpaOTsc.spatial_sc(sc_data = df_sc, is_data = df_is, is_dmat = is_dmat)
    issc.transport_plan(C**2, alpha = 0, rho = 1.0, epsilon = 0.1, cor_matrix = mcc, scaling = False)
        
    print('spatial')
        
    gamma = issc.gamma_mapping
    for j in range(gamma.shape[1]):
        gamma[:,j] = gamma[:,j]/np.sum(gamma[:,j])
    X_pred = np.matmul(gamma.T, np.array(issc.sc_data.values))
        
    result = pd.DataFrame(data = X_pred, columns = issc.sc_data.columns.values)
    test_genes = test_list
    result = result.loc[:, test_genes]
    print(result)
    return result

def stPlus_impute(K):
    global RNA_data, Spatial_data, outdir, train_gene, predict_gene
    
    test_list = predict_gene[K]
    train_list = train_gene[K]
    
    save_path_prefix = join(outdir, 'process_file/stPlus-demo')
    if not os.path.exists(join(outdir, "process_file")):
        os.mkdir(join(outdir, "process_file"))
    stPlus_res = stPlus(Spatial_data[train_list], RNA_data.T, test_list, save_path_prefix)
    
    return stPlus_res


def Tangram_impute(K, annotate = None, modes = 'clusters', density = 'rna_count_based'): # modified from tutorial
    import torch
    from torch.nn.functional import softmax, cosine_similarity, sigmoid
    import tangram as tg
    print ('We run Tangram for this \n' + DataDir)
    global RNA_data_adata, Spatial_data_adata, locations, train_gene, predict_gene
    

    test_list = predict_gene[K]
    test_list = [x.lower() for x in test_list]
    train_list = train_gene[K]
    
    locations = np.array(locations[['x_centroid','y_centroid']])
    
    if annotate == None:
        RNA_data_adata_label = RNA_data_adata
        sc.pp.normalize_total(RNA_data_adata_label)
        sc.pp.log1p(RNA_data_adata_label)
        sc.pp.highly_variable_genes(RNA_data_adata_label)
        RNA_data_adata_label = RNA_data_adata_label[:, RNA_data_adata_label.var.highly_variable]
        sc.pp.scale(RNA_data_adata_label, max_value=10)
        sc.tl.pca(RNA_data_adata_label)
        sc.pp.neighbors(RNA_data_adata_label)
        sc.tl.leiden(RNA_data_adata_label, resolution = 0.5)
        RNA_data_adata.obs['leiden']  = RNA_data_adata_label.obs.leiden
        tg.pp_adatas(RNA_data_adata, Spatial_data_adata, genes=train_list)
    else:
        global CellTypeAnnotate 
        RNA_data_adata.obs['leiden']  = CellTypeAnnotate
        tg.pp_adatas(RNA_data_adata, Spatial_data_adata, genes=train_list)
        
    device = torch.device('cuda:0')
    
    if modes == 'clusters':
        ad_map = tg.map_cells_to_space(RNA_data_adata, Spatial_data_adata, device = device, mode = modes, cluster_label = 'leiden', density_prior = density)
        ad_ge = tg.project_genes(ad_map, RNA_data_adata, cluster_label = 'leiden')
    else:
        ad_map = tg.map_cells_to_space(RNA_data_adata, Spatial_data_adata, device = device)
        ad_ge = tg.project_genes(ad_map, RNA_data_adata)
        
    test_list = list(set(ad_ge.var_names) & set(test_list))
    test_list = np.array(test_list)
    pre_gene = pd.DataFrame(ad_ge[:,test_list].X, index=ad_ge[:,test_list].obs_names, columns=ad_ge[:,test_list].var_names)
    
    return pre_gene


# # Data Input

# # First, you can download Example datasets from our websites 

# In[50]:


# test_cells = sys.argv[1]
# iterations = sys.argv[2]
# fraction = sys.argv[3]
# methods = list(sys.argv[4:])



methods = ['spage']
# methods = ['spage', 'novosparc', 'gimvi', 'tangram','spaotsc']
# methods = ['spaotsc', 'novosparc', 'gimvi', 'tangram']



test_cells = 10
iterations = 2
fraction = 10

methods = list(sys.argv[3:])
iterations = sys.argv[2]
test_cells = sys.argv[1]

print(iterations, methods)


# fraction = int(fraction)
iterations = int(iterations)
test_cells = int(test_cells)


# In[51]:


"""
@author: wen zhang
This function integrates two single-cell datasets, spatial and scRNA-seq, 
and predictes the expression of the spatially unmeasured genes from the scRNA-seq data.

Parameters
-------
RNA_file : str
    scRNA-seq data count file with Tab-delimited (cells X genes).
Spatial_file : str
    spatial count data file with Tab-delimited, please note that the file has no index.
location_file : str
    spatial spot coordinate file name with Tab-delimited, please note that the file has no index.
device : str
    Option,  ['CPU','GPU'], defaults to 'CPU'
train_gene : list
    genes for integrations, you can support more than one train list.
predict_gene : list
    genes for prediction, you can support more than one test list.
outdir : str
    result file stored direction    
"""

DataDir = 'xenium data/'

RNA_file = DataDir + 'Yao_150kcells_subsample_with_annotations.h5ad'

if fraction == 10:
    
    Spatial_file = DataDir + 'split_10/' + 'spatial_data.csv'
    location_file = DataDir + 'split_10/' + 'locations.csv'
    
else: 
    
    Spatial_file = DataDir + 'spatial_data.csv'
    location_file = DataDir + 'locations.csv'


device = 'GPU'



# In[ ]:





# In[52]:


outdir = 'xenium data/Figure2'

if not os.path.exists(outdir):
    os.mkdir(outdir)

outdir = outdir + '/downsampling_' + str(test_cells) 

if not os.path.exists(outdir):
    os.mkdir(outdir)
    
outdir = outdir + '/iterations_' + str(iterations)

if not os.path.exists(outdir):
    os.mkdir(outdir)
    
outdir = outdir + '/split_' + str(fraction)

if not os.path.exists(outdir):
    os.mkdir(outdir)
    
outdir = outdir + '/impute'

if not os.path.exists(outdir):
    os.mkdir(outdir)


# In[53]:


# test_cells = 0.01

test_cells = 1/ int(test_cells)
print(test_cells)

rna_cells = int(149955 * test_cells)
spatial_spots = int(82941 *test_cells)

# rna_cells=1000
# spatial_spots=1000


# ### RNA data

# In[54]:


RNA_data_adata = sc.read_h5ad(RNA_file)
RNA_data_adata


# In[55]:


RNA_data = pd.DataFrame(RNA_data_adata.X, columns = list(RNA_data_adata.var.index), index = list(RNA_data_adata.obs['sample_name']))
RNA_data


# In[56]:



# shuffle the DataFrame rows
RNA_data = RNA_data.sample(frac = 1)

RNA_data = RNA_data.T
RNA_data = RNA_data.iloc[:,:rna_cells] # downsampling

RNA_data


# In[57]:


# check for low-quality cells
# RNA_data[RNA_data.columns[(RNA_data != 0).sum(axis=0) > 2000]]


# In[58]:


RNA_data_adata = sc.AnnData(RNA_data.T)
RNA_data_adata


# ### Spatial data

# In[59]:


Spatial_data = pd.read_csv(Spatial_file, index_col=0)
Spatial_data


# In[60]:


# shuffle the DataFrame rows

if fraction != 10:
    Spatial_data = Spatial_data.sample(frac = 1)#.reset_index(drop=True)

    Spatial_data = Spatial_data.iloc[:spatial_spots,:] # downsampling
    Spatial_data


# In[61]:


Spatial_data.index


# In[62]:


# save to calculate metrics

# Spatial_data.to_csv(outdir + '/spatial_data.csv', header = 1, index = 1)
# print(outdir + '/spatial_data.csv')


# In[63]:


Spatial_data_adata = sc.AnnData(Spatial_data)
Spatial_data_adata


# In[64]:


locations = pd.read_csv(location_file, index_col=0)
locations = locations[locations['cell_code'].isin(list(Spatial_data.index))]
locations


# ### Training and testing list

# In[65]:


import math
import random 

if fraction == 10:
      
    train_gene = pd.read_csv(DataDir + 'split_10/' + 'train_list.csv', index_col=0).values.tolist()
    predict_gene = pd.read_csv(DataDir + 'split_10/' + 'test_list.csv', index_col=0).values.tolist()
    
    print('train list shape: ', np.array(train_gene).shape)
    print('test list shape: ', np.array(predict_gene).shape)
    
    train_gene_df = pd.DataFrame(train_gene)
    test_gene_df = pd.DataFrame(predict_gene)

    
else: 
    RNA_data_train = RNA_data

    RNA_data_train = RNA_data_train.loc[(RNA_data_train.sum(axis=1) != 0)] # avoiding error in SpaGE
    RNA_data_train = RNA_data_train.loc[(RNA_data_train.var(axis=1) != 0)]

    Spatial_data_cross = list(RNA_data_train.index&Spatial_data.columns) #list(Spatial_data.columns)

    n_total = len(Spatial_data_cross)
    print('total genes: ', n_total)

    n_test = int(len(Spatial_data_cross)/fraction) 
    # n_test = math.ceil(len(Spatial_data_cross)/fraction) 

    print('test genes: ', n_test)

    train_gene = []
    predict_gene = []

    # # random shuffle
    # for i in range(iterations):
    #     print(i)

    #     random.shuffle(Spatial_data_cross)

    #     Spatial_data_test = Spatial_data_cross[:n_test]
    #     Spatial_data_train = Spatial_data_cross[n_test:]

    #     train_gene.append(Spatial_data_train)
    #     predict_gene.append(Spatial_data_test)


    for i in range(iterations):
        print(i)

        Spatial_data_test = Spatial_data_cross[n_test*i : n_test*(i+1)]
        Spatial_data_train = [item for item in Spatial_data_cross if item not in Spatial_data_test]

        train_gene.append(Spatial_data_train)
        predict_gene.append(Spatial_data_test)


    print('train list shape: ', np.array(train_gene).shape)
    print('test list shape: ', np.array(predict_gene).shape)

     # save lists
    train_gene_df = pd.DataFrame(train_gene)
    test_gene_df = pd.DataFrame(predict_gene)

    # train_gene_df.to_csv(outdir + '/train_list.csv')
    # test_gene_df.to_csv(outdir + '/test_list.csv')


# # Predicting undetected transcripts by each method

# ### SpaGE

# In[66]:


list(range(len(train_gene)))[:iterations]


# In[69]:


if 'spage' in methods:
    
    savedir = outdir + '/SpaGE'

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    print(savedir)

    train_gene_df.to_csv(savedir + '/train_list.csv')
    test_gene_df.to_csv(savedir + '/test_list.csv')

    Spatial_data.to_csv(savedir + '/spatial_data.csv', header = 1, index = 1)
    locations.to_csv(savedir + '/locations.csv')
    
    KFOLD = list(range(len(train_gene)))[:iterations]
    
    print(KFOLD)
    
    for K in KFOLD:
    
        with multiprocessing.Pool(10) as pool:
            result = pd.concat(pool.map(SpaGE_impute, [K]), axis = 1)

            print(result)

        result.to_csv(savedir + '/SpaGE_impute'  + '_iteration' + str(K) + '.csv', header = 1, index = 1)

        print('saved file: ', savedir + '/SpaGE_impute' + '_iteration' + str(K)  + '.csv')


# ### SpaOTsc

# In[ ]:





# ### novoSpaRc

# In[ ]:





# # GPU Platform gimVI, Tangram, and stPlus

# ### gimVI

# In[21]:


if 'gimvi' in methods:

    savedir = outdir + '/gimVI'

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    train_gene_df.to_csv(savedir + '/train_list.csv')
    test_gene_df.to_csv(savedir + '/test_list.csv')

    Spatial_data.to_csv(savedir + '/spatial_data.csv', header = 1, index = 1)
    locations.to_csv(savedir + '/locations.csv')

    
    KFOLD = list(range(len(train_gene)))[:iterations]
    
    for K in KFOLD: 
        with multiprocessing.Pool(10) as pool:
            result = pd.concat(pool.map(gimVI_impute, [K]), axis = 1)
        result.to_csv(savedir + '/gimVI_impute' + '_iteration' + str(K) + '.csv', header = 1, index = 1)


# ### Tangram

# In[45]:


if 'tangram' in methods:
    
    savedir = outdir + '/Tangram'

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    train_gene_df.to_csv(savedir + '/train_list.csv')
    test_gene_df.to_csv(savedir + '/test_list.csv')

    Spatial_data.to_csv(savedir + '/spatial_data.csv', header = 1, index = 1)
    locations.to_csv(savedir + '/locations.csv')

    
    KFOLD = list(range(len(train_gene)))[:iterations]
    
    for K in KFOLD:
        with multiprocessing.Pool(10) as pool:
            result = pd.concat(pool.map(Tangram_impute, KFOLD), axis = 1)
        result.to_csv(savedir + '/Tangram_impute' + '_iteration' + str(K)  + '.csv', header = 1, index = 1)


# ## CPU slow

# ### novosparc

# In[46]:


if 'novosparc' in methods:
    savedir = outdir + '/novoSpaRc'

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    train_gene_df.to_csv(savedir + '/train_list.csv')
    test_gene_df.to_csv(savedir + '/test_list.csv')

    Spatial_data.to_csv(savedir + '/spatial_data.csv', header = 1, index = 1)
    locations.to_csv(savedir + '/locations.csv')

    
    KFOLD = list(range(len(train_gene)))[:iterations]
    
    for K in KFOLD:
        with multiprocessing.Pool(10) as pool:
            result = pd.concat(pool.map(novoSpaRc_impute, KFOLD), axis = 1)
        result.to_csv(savedir + '/novoSpaRc_impute' + '_iteration' + str(K)  + '.csv', header = 1, index = 1)
        print('saved file: ', savedir + '/novoSpaRC_impute' + '_iteration' + str(K)  + '.csv')


# ### Spaotsc

# In[ ]:


if 'spaotsc' in methods:
    
    savedir = outdir + '/SpaOTsc'

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    train_gene_df.to_csv(savedir + '/train_list.csv')
    test_gene_df.to_csv(savedir + '/test_list.csv')

    Spatial_data.to_csv(savedir + '/spatial_data.csv', header = 1, index = 1)
    locations.to_csv(savedir + '/locations.csv')

    
    KFOLD = list(range(len(train_gene)))[:iterations]
    print(KFOLD)
    
    for K in KFOLD:
        
        with multiprocessing.Pool(10) as pool:
            result = pd.concat(pool.map(SpaOTsc_impute, KFOLD), axis = 1)

        result.to_csv(savedir + '/SpaOTsc_impute' + '_iteration' + str(K)  + '.csv', header = 1, index = 1)
        print('saved file: ', savedir + '/SpaOTsc_impute' + '_iteration' + str(K)  + '.csv')



# In[ ]:





# In[ ]:





# ### stPlus

# In[ ]:


if 'stPlus' in methods:

    from stPlus import *
    # KFOLD = list(range(len(train_gene)))
    result = pd.DataFrame()
    for n in KFOLD:
        tmp = stPlus_impute(n)
        result = pd.concat([result, tmp], axis=1)
    result.to_csv(outdir +  '/stPlus_impute'  + '.csv',header = 1, index = 1)


# # Seurat impute file process

# In[18]:


# input_path = 'DataUpload/Dataset4/'
# output_path = 'FigureData/Figure2/Dataset4/'
# os.system('Rscript SpatialBenchmarking/FigureData/Figure2/RCodes/Seurat.r ' + input_path + ' ' + output_path)

# PATH = 'FigureData/Figure2/'
# DataSets = ['Dataset4']
# for Data in DataSets:
#     impute_count_file = PATH + Data + '/Seurat_impute.csv'
#     if os.path.exists(impute_count_file):
#         df = pd.read_csv(impute_count_file, header = 0, index_col = 0)
#         if df.index[0] == 0:
#             A = df
#             continue
#         print (impute_count_file)
#         Index = [(int(x)-1) for x in df.index]
#         Index = sorted(Index)
#         df.index = Index
#         B = df
#         df.to_csv(impute_count_file)
#     else:
#         print ('This outdir is none for Seurat : ' + PATH + Data)


# # LIGER impute file process

# In[ ]:


# input_path = 'DataUpload/Dataset4/'
# output_path = 'SpatialBenchmarking/FigureData/Figure2/Dataset4/'

# os.system('Rscript SpatialBenchmarking/FigureData/Figure2/RCodes/LIGER.r'  + 'DataUpload/Dataset4/' + ' ' + 'SpatialBenchmarking/FigureData/Figure2/Dataset4/')



# In[23]:


# input_path = 'DataUpload/Dataset4/'
# output_path = 'SpatialBenchmarking/FigureData/Figure2/Dataset4/'
# os.system('Rscript SpatialBenchmarking/FigureData/Figure2/RCodes/LIGER.r ' + input_path + ' ' + output_path)

# PATH = 'SpatialBenchmarking/FigureData/Figure2/Dataset4/'
# DataSets = ['Dataset4']
# for Data in DataSets:
#     impute_count_file = PATH + Data + '/LIGER_impute.csv'
#     if os.path.exists(impute_count_file):
#         df = pd.read_csv(impute_count_file, header = 0, index_col = 0)
#         if df.index[0] == 0:
#             A = df
#             continue
#         print (impute_count_file)
#         Index = [(int(x.replace('V', '')) -1) for x in df.index]
#         Index = sorted(Index)
#         df.index = Index
#         B = df
#         df.to_csv(impute_count_file)
#     else:
#         print ('This outdir is none for LIGER : ' + PATH + Data)
        


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




