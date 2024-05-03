import warnings
warnings.filterwarnings("ignore") 
import random
import scipy.sparse as sparse
from scipy.sparse import csr_matrix, issparse
from Banksy_py.banksy.initialize_banksy import initialize_banksy
from Banksy_py.banksy.run_banksy import run_banksy_multiparam



## modify the X coords of each sample so that they are processed independently
def adapt_banksy_for_multisample(adata):
    adata.obs['Xres']=adata.obs['x_centroid']
    adata.obs['Yres']=adata.obs['y_centroid']
    gap=adata.obs['x_centroid'].max()+int(adata.obs['x_centroid'].max()/10)
    samplekey='sample'
    counter=0
    sampsel=[]
    for s in adata.obs[samplekey].unique():
        sampsel.append(s)
        adata.obs.loc[~adata.obs[samplekey].isin(sampsel),'Xres']=adata.obs.loc[~adata.obs[samplekey].isin(sampsel),'Xres']+gap
        counter=counter+1
    adata.obs['Yres']=adata.obs['Yres']+1
    adata.obsm['spatial_sample_specific']=adata.obsm['spatial']
    adata.obsm['spatial']=np.array(adata.obs.loc[:,['Xres','Yres']])
    return adata
def define_palette(n_colors=50):
    from random import randint
    colorlist = []
    n = n_colors
    for i in range(n):
        colorlist.append('#%06X' % randint(0, 0xFFFFFF))
    return colorlist

def domains_by_banksy(adata,plot_path:str,banksy_params:dict):
    adata=adapt_banksy_for_multisample(adata)
    coord_keys = ('Xres', 'Yres','spatial')
    prev_clusters=[e for e in adata.obs.columns if clustering_params['clustering_alg'] in e]
    if len(prev_clusters)>0:
        annotation_key=prev_clusters[0]
    else:
        adata.obs['default_clustering']='0'
        adata.obs['default_clustering']=adata.obs['default_clustering'].astype('category')
        annotation_key='default_clustering'

    banksy_dict = initialize_banksy(
        adata,
        coord_keys,
        banksy_params['k_geom'],
        nbr_weight_decay=banksy_params['nbr_weight_decay'],
        max_m=banksy_params['max_m'],
        plt_edge_hist=True,
        plt_nbr_weights=True,
        plt_agf_angles=False,
        plt_theta=False)
    results_df = run_banksy_multiparam(
        adata,
        banksy_dict,
        banksy_params['lambda_list'], banksy_params['resolutions'],
        color_list = define_palette(n_colors=100), max_m = max_m, filepath = plot_path,
        key = coord_keys, pca_dims = banksy_params['pca_dims'],
        annotation_key = annotation_key, max_labels = None,
        cluster_algorithm = banksy_params['cluster_algorithm'], match_labels = False, savefig = save, add_nonspatial = False,
        variance_balance = False,
    )
    adata_res=results_df.loc[results_df.index[0],'adata']
    res=adata_res.obs.loc[:,['unique_cell_id',results_df.index[0]]]
    res.columns=['unique_cell_id','spatial_domain']
    adata_res=results_df.loc[results_df.index[0],'adata']
    res=adata_res.obs.loc[:,['unique_cell_id',results_df.index[0]]]
    res.columns=['unique_cell_id','spatial_domain']
    id2domain=dict(zip(res['unique_cell_id'],res['spatial_domain']))
    adata.obs['banksy_domain']=adata.obs['unique_cell_id'].map(id2domain).astype(str)
    adata.obsm['spatial']=adata.obsm['spatial_sample_specific']
    return adata,adata_res
    
    
def format_data_neighs_colapse(adata,condit,neighs=10):
    try:
        adata.obsm['spatial']
    except:
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

def domains_by_rbd(adata,hyperparameters_rbd:dict): #read-based domains
    sample_key='sample'  # sample is considered the "sample_key"
    anndata_list=[]
    for sample in adata.obs[sample_key].unique():
        adata_copy_int = adata[adata.obs[sample_key]==sample]
        adataneigh2=format_data_neighs_colapse(adata,sample_key,neighs=hyperparameters_rbd['neighbors'])
        anndata_list.append(adataneigh2)
    adataneigh=sc.concat(anndata_list) 
    adataneigh.X=np.nan_to_num(adataneigh.X)
    adataneigh=adataneigh[adataneigh.obs['total_counts']>3]
    adataneigh.raw=adataneigh
    sc.pp.neighbors(adataneigh, n_neighbors=20,n_pcs=0)
    sc.tl.umap(adataneigh,min_dist=0.1)
    if hyperparameters_rbd['clustering_algorithm']=='leiden':
        sc.tl.leiden(adataneigh,resolution=hyperparameters_rbd['resolution'],key_added='rbd_domain')
    if hyperparameters_rbd['clustering_algorithm']=='louvain':
        sc.tl.louvain(adataneigh,resolution=hyperparameters_rbd['resolution'],key_added='rbd_domain')    
    if hyperparameters_rbd['clustering_algorithm']=='leiden':
        sc.tl.leiden(adataneigh,resolution=hyperparameters_rbd['resolution'],key_added='rbd_domain')
    if hyperparameters_rbd['clustering_algorithm']=='louvain':
        sc.tl.louvain(adataneigh,resolution=hyperparameters_rbd['resolution'],key_added='rbd_domain')    
    id2domain=dict(zip(adataneigh.obs['unique_cell_id'],adataneigh.obs['rbd_domain']))
    adata.obs['rbd_domain']=adata.obs['unique_cell_id'].map(id2domain).astype(str)
    return adata,adataneigh


def format_data_neighs(adata,sname,condit,neighs=10):
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
    adata1.obs['sample']=condit
    adata1.obs['condition']=condit
    return adata1
    
def compare_domains(adata,domain_keys:list,save=True,plot_path='./'):
    import sklearn.metrics as sk
    ARI=pd.DataFrame(index=domain_keys,columns=domain_keys)
    for ind in ARI.index:
        for col in ARI.columns:
            ARI.loc[ind,col]=sk.adjusted_rand_score(adata.obs[ind],adata.obs[col])
    ARI.index=[i.replace('_domain','') for i in ARI.index]
    ARI.columns=[i.replace('_domain','') for i in ARI.columns]
    plt.figure()
    sns.clustermap(ARI.astype(float),cmap='Greens',figsize=(len(ARI.columns)*1.5,len(ARI.columns)*1.5),dendrogram_ratio=0.15)
    if save==True:
        plt.savefig(plot_path+'clustermap_ARI_domains.pdf')
    for s in adata.obs['sample'].unique():
            fig,ax=plt.subplots(nrows=len(domain_keys))
            adatasub=adata[adata.obs['sample']==s]
            s=0
            for groupby in domain_keys:
                sc.pl.spatial(adatasub,color=groupby,spot_size=40,ax=ax[s],show=False)
                s=s+1
            if save==True:
                plt.savefig(plot_path+'map_all_domains_'+str(s)+'.pdf')
    
    return ARI