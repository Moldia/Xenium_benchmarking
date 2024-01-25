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

def missegmentation_simulation(adata_sc_sub,missegmentation_percentage=0.1):
    exp=adata_sc_sub.to_df()
    if missegmentation_percentage>0:
        cells_affected=int(exp.shape[0]*(missegmentation_percentage/100))
        for num in tqdm(range(0,cells_affected)):
            missegmentation_importance=random.sample(range(0,10),1)[0]*0.1
            source=random.sample(list(exp.index),1)[0]
            target=random.sample(list(exp.index),1)[0]
            exp.loc[target,:]=exp.loc[target,:]+(exp.loc[source,:]*missegmentation_importance)
        adata_sc_sub.X=np.array(exp.astype(int))
    return adata_sc_sub
def noise_adder(adata_sc,percentage_of_noise=0.1):
    noise_events=int(np.sum(adata_sc.X)*(percentage_of_noise/100))
    for i in range(0,noise_events):
        x=random.sample(range(0,adata_sc.X.shape[0]),1)
        y=random.sample(range(0,adata_sc.X.shape[1]),1)
        change=random.sample([-1,1],1)
        adata_sc.X[x,y]=adata_sc.X[x,y]+change
    adata_sc.X[adata_sc.X<-1]=0
    return adata_sc

def subset_of_single_cell(adata_sc_sub,markers,random_markers_percentage=0,
                          reads_x_cell=None,number_of_markers=200,
                         n_reads_x_gene=40,percentage_of_noise=0.1,ms_percentage=0.1):
    mk=[]
    number_of_markers_x_cluster=int(np.ceil(number_of_markers/markers.shape[1])+1)
    for ind in  markers.index[0:number_of_markers_x_cluster]:
        for col in markers.columns:
            mk.append(markers.loc[ind,col])
    mk=np.unique(mk)
    if len(mk)>number_of_markers:
        mk=random.sample(mk,number_of_markers)
    
    sampling_amount=(int(len(np.unique(mk))*(random_markers_percentage/100)))
    newg=random.sample(list(adata_sc_sub.var.index),sampling_amount)
    sel=random.sample(range(0,len(mk)),len(newg))
    ii=0
    for el in sel:
        mk[el]=newg[ii]
        ii=ii+1
    mk=list(mk)
    adata_sc=adata_sc_sub[:,adata_sc_sub.var.index.astype(int).isin(mk)].copy()
     ###THIS STEP MAKES MORE EXPRESSED GENES TO BE LESS DETECTED 
    norm_factors=np.percentile(adata_sc.X,99.9,axis=0)+1
    
    #number of reads per gene (aprox)
   ##WE HAVE COMENTED THIS# adata_sc.X=np.array(adata_sc.to_df().div(np.array(norm_factors)/n_reads_x_gene,axis=1).astype(int))
    # segmentation_simulator
    adata_sc=missegmentation_simulation(adata_sc,missegmentation_percentage=ms_percentage)
    ###ADD RANDOM NOISE
    adata_sc=noise_adder(adata_sc,percentage_of_noise=percentage_of_noise)
    if reads_x_cell is not None:
        sc.pp.downsample_counts(adata_sc,counts_per_cell=reads_x_cell)
    print(adata_sc.shape)
    return adata_sc
from sklearn.metrics import mutual_info_score

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

from sklearn.metrics import fowlkes_mallows_score

def compute_fmi(ground_truth, predicted):
    fmi_score = fowlkes_mallows_score(ground_truth, predicted)
    return fmi_score
def keep_nuclei_and_quality(adata1,overlaps_nucleus=1,qvmin=20):
    if overlaps_nucleus==1:
        subset1=adata1.uns['spots'].loc[adata1.uns['spots']['overlaps_nucleus']==overlaps_nucleus,:]
    if overlaps_nucleus==0:
        subset1=adata1.uns['spots']
    subset1=subset1[subset1['qv']>qvmin]
    ct1=pd.crosstab(subset1['cell_id'],subset1['feature_name'])
    adataobs=adata1.obs.loc[adata1.obs['cell_id'].isin(ct1.index),:]
    av=adata1.var.index[adata1.var.index.isin(ct1.columns)]#.isin(adata1.var.index)
    ct1=ct1.loc[:,av]
    adataobs.index=adataobs['cell_id']
    adataobs.index.name='ind'
    ct1=ct1.loc[ct1.index.isin(adataobs['cell_id']),:]
    adata1nuc=sc.AnnData(np.array(ct1),obs=adataobs)#,var=adata1.var)
    return adata1nuc

# using a baseline, modify parameters one by one
# making it into a function
def allcombs(adata):
    default={'nuc':1,'npcs':0,'neighs':16,'qvs':0,'scale':True,'hvgs':False,'ts':100,'norm':True,'lg':True} 
    try:
        del allres
    except:
        print('first run')
    nucdic={1:'nuc',0:'cell'}
    neighs=[16,12,6,20]
    npcs=[0,15,25]
    ts=[100,10,1000,None]
    qvs=[0,20,30]
    scales=[True,False]
    hvgs=[False,True]
    norm=[True,False]
    logs=[True,False]
    perc=0.1
    for qv in qvs:
        for nuc in [1,0]:
            if nuc==1:
                adata_f1=keep_nuclei_and_quality(adata.copy(),overlaps_nucleus=1,qvmin=qv)
                adata_f1=adata_f1.copy()
            if nuc==0:
                adata_f1=keep_nuclei_and_quality(adata.copy(),overlaps_nucleus=0,qvmin=qv)
                adata_f1.uns['spots']=[]
            for target_sum in ts:
                for neigh in neighs:
                    for nm in norm:
                        for npc in npcs:
                            for scale in scales:
                                for lgo in logs:
                                    for hvg in hvgs:
                                        nv=neigh!=default['neighs']
                                        nv=nv+(nuc!=default['nuc'])
                                        nv=nv+(qv!=default['qvs'])
                                        nv=nv+(npc!=default['npcs'])
                                        nv=nv+(scale!=default['scale'])
                                        nv=nv+(hvg!=default['hvgs'])
                                        nv=nv+(target_sum!=default['ts'])
                                        nv=nv+(nm!=default['norm'])
                                        nv=nv+(lgo!=default['lg'])
                                        #print(nv)
                                        if nv==0:
                                            adata_f2=adata_f1.copy()
                                            print('DEFAULT')
                                            result,ttc=main_preprocessing(adata_f2,neigh=neigh,npc=npc,target_sum=target_sum,scale=scale,hvg=hvg,default=True,norm=nm,lg=lgo)
                                            sumresult=result.obs.loc[:,['leiden_1_4','louvain_1_4']]
                                            sumresult.columns=['DEFAULT_leid','DEFAULT_louv']
                                            try:
                                                allres=pd.merge(allres,sumresult, left_index=True, right_index=True)
                                            except:
                                                allres=sumresult
                                        if nv==1:
                                            adata_f2=adata_f1.copy()
                                            print(nucdic[nuc]+'_qv'+str(qv)+'_neighs'+str(neigh)+'_npc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale))
                                            result=main_preprocessing(adata_f2,neigh=neigh,npc=npc,target_sum=target_sum,scale=scale,hvg=hvg,total_clusters=ttc,norm=nm,lg=lgo)
                                            sumresult=result.obs.loc[:,['leiden_1_4','louvain_1_4']]
                                            sumresult.columns=[nucdic[nuc]+'_qv'+str(qv)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale)+'_norm_'+str(nm)+'_log_'+str(lgo)+'_leid',nucdic[nuc]+'_qv'+str(qv)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale)+'_norm_'+str(nm)+'_log_'+str(lgo)+'_louv']
                                            try:
                                                allres=pd.merge(allres,sumresult, left_index=True, right_index=True)
                                            except:
                                                allres=sumresult
    return allres



def main_preprocessing(adata,target_sum=100,mincounts=10,mingenes=3,neigh=15,npc=0,nuc=1,scale=False,hvg=False,default=False,total_clusters=30,default_resol=1.6,logstatus=True,normstatus=True):
#    print(adata.shape)
    adata.layers['raw']=adata.X.copy()
    sc.pp.filter_cells(adata,min_counts=mincounts)
    sc.pp.filter_cells(adata,min_genes=mingenes)
    adata.raw=adata
    adata.layers['raw']=adata.X.copy()
    if normstatus=='pearson':
        sc.experimental.pp.normalize_pearson_residuals(adata)
    else:
        if normstatus==True:
            sc.pp.normalize_total(adata, target_sum=target_sum)
        if logstatus==True:
            sc.pp.log1p(adata)
        if hvg==True:
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=12, min_disp=0.25)
            adata = adata[:, adata.var.highly_variable]
        if scale==True:
            sc.pp.scale(adata)
    #remove nans, inf and -ind
    adata.X=np.array(adata.to_df().fillna(0))
    adata.X[adata.X==-np.inf]=0
    adata.X[adata.X==np.inf]=0
#    adata.X.fillna(0, inplace=True)
    sc.pp.pca(adata)
    print(adata.shape)
    sc.pp.neighbors(adata, n_neighbors=neigh, n_pcs=npc)
#    sc.tl.leiden(adata,resolution=2.2,key_added='leiden_2_2')
    resol=default_resol
    sc.tl.leiden(adata,resolution=resol,key_added='leiden_1_4')
    numclust=int(np.max(adata.obs['leiden_1_4'].astype(int)))
    targetnum=total_clusters
    if default==True:
        targetnum=numclust
    if abs(numclust-targetnum)>3:
        i=0
        while abs(int(numclust)-int(targetnum))>3:
       #     print('iter '+str(i))
            if numclust>targetnum:
                resol=resol-0.05
            if numclust<targetnum:  
                resol=resol+0.05
            sc.tl.leiden(adata,resolution=resol,key_added='leiden_1_4')
            numclust=np.max(adata.obs['leiden_1_4'].astype(int))
            i=i+1
            if i >100:
                break
    resol=default_resol
    sc.tl.louvain(adata,resolution=resol,key_added='louvain_1_4')
    numclust=np.max(adata.obs['louvain_1_4'].astype(int))
    if abs(numclust-targetnum)>3:
        i=0
        while abs(int(numclust)-int(targetnum))>3:
    #        print('iter '+str(i))
            if numclust>targetnum:
                resol=resol-0.05
            if numclust<targetnum:  
                resol=resol+0.05
            sc.tl.leiden(adata,resolution=resol,key_added='louvain_1_4')
            numclust=np.max(adata.obs['louvain_1_4'].astype(int))
            i=i+1
            if i >100:
                break
#    sc.tl.umap(adata)
#    sc.pl.umap(adata,color=['leiden_1_4','louvain_1_4','leiden'],ncols=1)
    if default==True:
        return adata,targetnum
    else:
        return adata
    
# using a baseline, modify parameters one by one
# making it into a function
################################################DEFAULT KEY################################
def allcombs_simulated(adata,default_key='class'):
    default={'nuc':1,'npcs':0,'neighs':12,'qvs':0,'scale':False,'hvgs':False,'ts':1000000}# default is never done
    try:
        del allres
    except:
        print('first run')
    tot=0
    total_combs=3000#3000 
    neighs=[12,16]#12
    npcs=[0,20,40]
    ts=[10,100,1000,None]#1000
    qvs=[0]
    scales=[False,True]
    hvgs=[False,True]#True
    log=[True,False]#100,80
    norm=[True,False,'pearson']#,
    perc=0.1
    df_resol=0.4
    nuc=1
    for lg in log:
        for target_sum in ts:
            for neigh in neighs:
                for npc in npcs:
                    for scale in scales:
                        for hvg in hvgs:
                            for nm in norm:
                                if tot>total_combs:
                                    print('Broken')
                                    break
                                tot=tot+1
                                print(tot)
                                try:
                                    if nm=='pearson':
                                        print(str(lg)+str(scale)+str(hvg)+str(ts))
                                        if lg==True and scale==True and hvg==True and target_sum==None:
                                            ttc=len(np.unique(adata.obs[default_key]))
                                            adata_f1=adata.copy()
                                            adata_f2=adata_f1.copy()
                                            result=main_preprocessing(adata_f2,neigh=neigh,npc=npc,target_sum=target_sum,scale=scale,hvg=hvg,total_clusters=ttc,default_resol=df_resol,logstatus=lg,normstatus=nm)
                                            sumresult=result.obs.loc[:,['leiden_1_4','louvain_1_4']]
                                            sumresult.columns=['norm'+str(nm)+'_lg'+str(lg)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale)+'_leid','norm'+str(nm)+'_lg'+str(lg)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale)+'_louvain']
                                            try:
                                                allres=pd.merge(allres,sumresult, left_index=True, right_index=True)
                                            except:
                                                allres=sumresult
                                    else:
                                        ttc=len(np.unique(adata.obs[default_key]))
                                        adata_f1=adata.copy()
                                        adata_f2=adata_f1.copy()
                                        result=main_preprocessing(adata_f2,neigh=neigh,npc=npc,target_sum=target_sum,scale=scale,hvg=hvg,total_clusters=ttc,default_resol=df_resol,logstatus=lg,normstatus=nm)
                                        sumresult=result.obs.loc[:,['leiden_1_4','louvain_1_4']]
                                        sumresult.columns=['norm'+str(nm)+'_lg'+str(lg)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale)+'_leid','norm'+str(nm)+'_lg'+str(lg)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale)+'_louvain']
                                        try:
                                            allres=pd.merge(allres,sumresult, left_index=True, right_index=True)
                                        except:
                                            allres=sumresult
                                except:
                                    print('norm'+str(nm)+'_lg'+str(lg)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale)+'_leid','lg'+str(lg)+'_ng'+str(neigh)+'_pc'+str(npc)+'_ts'+str(target_sum)+'_hvg_'+str(hvg)+'_scale_'+str(scale))

    allres['RESULTS']=adata.obs[default_key]                            
    return allres
