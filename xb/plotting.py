import numpy as np
import anndata as ad
from anndata import AnnData
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import random

def map_of_clusters(adata,key='leiden',clusters='all',size=8,background='white',figuresize=(10,7),save=None,format='pdf'):
    
   """ Make spatial plots based on a given adata object.
   
    Parameters:
    key (str): the terms in adata.obs that you want to plot
    clusters (str or list):'all' for plotting all clusters in a single plot, 'individual': for plots of individual genes, or ['3','5'] (your groups between square brackets to plot only some clusters
    size: to change the size of your spots
    background (str): color of the background 
    figuresize (tupple): to specify the size of your figure
    save (boolean or str): whether want to save your figure. If so, please add the PATH of the folder where you want to save it
    format (str): specify the format in which you want to save your figure (i.e. '.pdf', '.png')

    Returns:
    None
   
   """



    
    try:
        adata.obs[key]=adata.obs[key].astype(int)
        colors=dict(zip(np.unique(adata.obs[key]),adata.uns[key+'_colors']))
    except:
        colors=dict(zip(np.unique(adata.obs[key]),adata.uns[key+'_colors']))
    #cl.apply(lambda x: colors[x])
    plt.rcParams['figure.facecolor'] = background
    if clusters=='all':
        cl=adata.obs[key]
        plt.figure(figsize=figuresize)
        figa=plt.scatter(x=adata.obs.X,y=adata.obs.Y,c=cl.apply(lambda x: colors[x]),s=size,linewidths=0, edgecolors=None)
        plt.axis('off')
        if not save==None:
            plt.savefig(save +'/map_all_clusters_'+str(size)+'_'+background+'_'+key+'.'+format)
    elif clusters=='individual':
        cl=adata.obs[key]
        for each in adata.obs[key].unique():
            adatasub=adata[adata.obs[key]==each]
            plt.figure(figsize=figuresize)
            plt.scatter(x=adata.obs.X,y=adata.obs.Y,c='grey',s=size/5,linewidths=0, edgecolors=None)
            cl=adatasub.obs[key]
            plt.scatter(x=adatasub.obs.X,y=adatasub.obs.Y,c=cl.apply(lambda x: colors[x]),s=size,linewidths=0, edgecolors=None)
            plt.axis('off')
            plt.title('Group: '+ str(each))
            if not save==None:
                plt.savefig(save +'/map_inidivdual_cluster_'+str(each)+'_'+str(size)+background+'_'+key+'.'+format)
    else:
        adatasub=adata[adata.obs[key].isin(clusters)]
        plt.figure(figsize=figuresize)
        plt.scatter(x=adata.obs.X,y=adata.obs.Y,c='grey',s=size/5,linewidths=0, edgecolors=None)
        cl=adatasub.obs[key]
        plt.scatter(x=adatasub.obs.X,y=adatasub.obs.Y,c=cl.apply(lambda x: colors[x]),s=size,linewidths=0, edgecolors=None)
        plt.axis('off')
        plt.legend()
        if not save==None:
                s=''
                for element in clusters:
                    s=s+str(element)
                print(s)
                plt.savefig(save +'/map_group_of_clusters_'+str(s)+'_'+str(size)+background+'_'+key+'.'+format)
#        plt.title('Group: '+ paste(clusters))

def generate_hex_colors(num_colors=70):
    """ Generate a list of hex colors.
    
    Parameters:
    num_colors(int): number of colors to generate
    
    Returns:
    hex_colors (list):list of randomly generated colors
   
   """   
    
    hex_colors = []
    for _ in range(num_colors):
        # Generate a random hex color
        color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
        hex_colors.append(color)
    return hex_colors



def plot_cell_counts(adata,plot_path:str,save=True,clustering_params={}):
      """ Plot the histogram of the counts detected per cell
   
    Parameters:
    adata (AnnData): AnnData object with the information of cells profiled
    plot_path (str): path where to save the generated plot, if needed
    save (boolean): whether to save or not the output path
    clustering_params (dict): list of parameters used for preprocessing and clustering the experiment 

    Returns:
    None
   
   """
    fig,ax=plt.subplots(ncols=2,figsize=(12,3),dpi=100)
    vivi=ax[0].hist(adata.obs['total_counts'],bins=200,color='#f9debd')
    ax[0].set_xlabel('Counts/cell')
    ax[0].set_ylabel('Total cells')
    ax[0].axvline(x = clustering_params['min_counts_x_cell'], color = '#786e8a', label = 'axvline - full height')
    vivi=ax[1].hist(adata.obs['expressed_genes'],bins=70,color='#f9debd')
    ax[1].set_xlabel('Genes/cell')
    ax[1].set_ylabel('Total cells')
    ax[1].axvline(x=clustering_params['min_genes_x_cell'], color = '#786e8a', label = 'axvline - full height')
    if save==True:
        plt.savefig(plot_path+'cell_counts_histogram.png',dpi=200)
        
        
def plot_domains(adata,groupby='nbd_domain'):
    """ Generate the spatial plots of the domains previously identified
   
    Parameters:
    adata (AnnData): AnnData object with the information of cells profiled
    groupby (str): Name of the column in adata.obs where the domain information is stored

    Returns:
    None

    """
    for s in adata.obs['sample'].unique():
        adatasub=adata[adata.obs['sample']==s]
        sc.pl.spatial(adatasub,color=groupby,spot_size=40)
