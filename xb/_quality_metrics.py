import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import ConvexHull, convex_hull_plot_2d


def cell_density(adata_sp: AnnData,pipeline_output=True):
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
    density=(adata_sp.shape[0]/hull.area)#*10e6
    return density 

def proportion_of_assigned_reads(adata_sp: AnnData,pipeline_output=True):
    """Proportion of assigned reads
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    proportion_assigned : float
       Proportion of reads assigned to cells / all reads decoded
    """   
    proportion_assigned=np.sum(adata_sp.layers['raw'])/adata_sp.uns['spots'].shape[0]
    return proportion_assigned


def median_reads_cells(adata_sp: AnnData,pipeline_output=True):
    """Median number of reads/cells
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    median_cells : float
       Median_number_of_reads_x_cell
    """   
    median_cells=np.median(np.sum(adata_sp.layers['raw'],axis=1))
    return median_cells

def mean_reads_cells(adata_sp: AnnData,pipeline_output=True):
    """Mean number of reads/cells
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    mean_cells : float
       Mean_number_of_reads_x_cell
    """   
    mean_cells=np.mean(np.sum(adata_sp.layers['raw'],axis=1))
    return mean_cells

def number_of_genes(adata_sp: AnnData,pipeline_output=True):
    """ Size of the gene panel present in the spatial dataset
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    number_of genes : float
       Number of genes present in the spatial dataset
    """   
    number_of_genes=adata_sp.shape[1]
    return number_of_genes

def number_of_cells(adata_sp: AnnData,pipeline_output=True):
    """ Number of cells present in the spatial dataset
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    number_of cells : float
       Number of cells present in the spatial dataset
    """   
    number_of_cells=adata_sp.shape[0]
    return number_of_cells

def percentile_5th_reads_cells(adata_sp: AnnData,pipeline_output=True):
    """5th percentile of number of reads/cells in the spatial experiment
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    median_cells : float
       Median_number_of_reads_x_cell
    """   
    pctile5=np.percentile(np.squeeze(np.asarray(np.sum(adata.layers['raw'],axis=1))),5)
    return pctile5

def mean_genes_cells(adata_sp: AnnData,pipeline_output=True):
    """Mean number of genes/cell in the spatial experiment
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    median_cells : float
       Mean number of genes per cell
    """   
    mean_genesxcell=np.mean(np.sum((adata_sp.layers['raw']>0)*1,axis=1))
    return mean_genesxcell

def percentile_95th_genes_cells(adata_sp: AnnData,pipeline_output=True):
    """Percentile 95 of genes/cell in the spatial experiment
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    median_cells : float
       Percentile 95 of genes per cell
    """   
    percentile95_genesxcell=np.percentile(np.sum((adata_sp.layers['raw']>0)*1,axis=1),95)
    return percentile95_genesxcell

def percentile_5th_genes_cells(adata_sp: AnnData,pipeline_output=True):
    """Percentile 5 of genes/cell in the spatial experiment
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    median_cells : float
       Percentile 5 of genes per cell
    """   
    percentile5_genesxcell=np.percentile(np.sum((adata_sp.layers['raw']>0)*1,axis=1),5)
    return percentile5_genesxcell

def median_genes_cells(adata_sp: AnnData,pipeline_output=True):
    """Median of genes/cell in the spatial experiment
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    median_cells : float
       Median of genes per cell
    """   
    median_genesxcell=np.median(np.sum((adata_sp.layers['raw']>0)*1,axis=1))
    return median_genesxcell




def percentile_95th_reads_cells(adata_sp: AnnData,pipeline_output=True):
    """5th percentile of number of reads/cells in the spatial experiment
    Parameters
    ----------
    adata_sp : AnnData
        annotated ``AnnData`` object with counts from spatial data
    pipeline_output : float, optional
        Boolean for whether to use the 
    Returns
    -------
    median_cells : float
       Median_number_of_reads_x_cell
    """   
    pctile95=np.percentile(np.squeeze(np.asarray(np.sum(adata.layers['raw'],axis=1))),95)
    return pctile95
    
