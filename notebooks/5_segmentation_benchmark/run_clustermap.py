#!/usr/bin/env python

import argparse
from collections import OrderedDict
from ClusterMap.clustermap import ClusterMap
import tifffile
import numpy as np
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign molecules to cells using ClusterMap')
    parser.add_argument('-m', '--molecules', required=True, type=str, 
        help='Input csv file in format [Gene, x, y]')
    parser.add_argument('-d', '--data', required=True, type=str, 
        help='Ouput data directory- should also contain segmented image')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input image file')
    parser.add_argument('-p', '--hyperparams', default=None, type=str,
        help='Dictionary of hyperparameters') 
    parser.add_argument('-id', '--id_code', required=True, type = str,
        help='ID of method to be used for saving')
    
    args = parser.parse_args()

    molecules = args.molecules
    data = args.data
    image = args.input
    hyperparams = eval(args.hyperparams)
    id_code = args.id_code

    #Make sure parameter dictionaries are not `None`
    if hyperparams is None: hyperparams = {}
    if hyperparams.get('model') is None: hyperparams['model'] = {'xy_radius':15}
    if hyperparams['model'].get('xy_radius') is None: hyperparams['model']['xy_radius']=15
    if hyperparams.get('preprocess') is None: hyperparams['preprocess'] = {}
    if hyperparams.get('segmentation') is None: hyperparams['segmentation'] = {}

    #Read and format input data
    dapi = tifffile.imread(image)
    num_dims=len(dapi.shape)
    spots = pd.read_csv(molecules)

    spots.rename(columns = {spots.columns[0]:'gene_name',
                        spots.columns[1]:'spot_location_1',
                        spots.columns[2]:'spot_location_2'
                        } , inplace = True)

    #Use gene id numbers instead of names
    genes, ids = np.unique(spots['gene_name'], return_inverse=True)
    spots['gene'] = ids+1
    spots = spots.astype({'spot_location_1':int, 'spot_location_2':int})
    gene_list=np.unique(ids)+1
    genes = pd.DataFrame(genes)

    #Create Model
    model = ClusterMap(spots=spots,dapi=dapi, gene_list=gene_list, num_dims=num_dims,z_radius=0,
        **(hyperparams['model']))

    #Preprocess
    model.preprocess(**(hyperparams['preprocess']))

    #Segment cells
    model.min_spot_per_cell = 2
    model.segmentation(**(hyperparams['segmentation']))

    #The original spots file is modified by Clustermap to include assignments
    #Copy and return spots-to-cell assignment
    assignments = spots.copy()    
    assignments.rename(columns = {'gene_name':'Gene',
                        'spot_location_1':'x',
                        'spot_location_2':'y',
                        'clustermap':'cell',
                        } , inplace = True)
    assignments.cell = assignments.cell+1

    #Save to csv
    assignments.to_csv(f'{data}/assignments_clustermap-{id_code}.csv', index = False)