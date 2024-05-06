from anndata import AnnData
import numpy as np
import pandas as pd
from ._quality_metrics import *

def all_quality_metrics(
    adata_sp: AnnData,
) -> pd.DataFrame:

    """ Compute all metrics quality metrics given an AnnData object 
    
    Parameters:
    adata_sp (AnnData): AnnData object with the cells of the experiment 

    Returns:
    resulting_metrics(DataFrame): Quality metrics computed for the dataset used as an input

   """

    
    #Generate metrics
    metrics = {}
    
    metrics['cellular_density']=cell_density(adata_sp)
    metrics['prop_reads_assigned']=proportion_of_assigned_reads(adata_sp)
    metrics['median_readsxcell']=median_reads_cells(adata_sp)
    metrics['mean_readsxcell']=mean_reads_cells(adata_sp)
    metrics['number_of_genes']=number_of_genes(adata_sp)
    metrics['number_of_cells']=number_of_cells(adata_sp)
    metrics['pct5_readsxcell']=percentile_5th_reads_cells(adata_sp)
    metrics['mean_genesxcell']=mean_genes_cells(adata_sp)
    metrics['pct95_genesxcell']=percentile_95th_genes_cells(adata_sp)
    metrics['pct5_genesxcell']=percentile_5th_genes_cells(adata_sp)
    metrics['median_genexcell']=median_genes_cells(adata_sp)
    metrics['pct95_readsxcell']=percentile_95th_reads_cells(adata_sp)
    
    resulting_metrics=pd.DataFrame.from_dict(metrics, orient='index')
    
    
    
    return resulting_metrics

