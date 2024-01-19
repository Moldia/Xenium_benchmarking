from typing import Union
from sklearn.cluster import MiniBatchKMeans as KMeans
from sklearn.preprocessing import normalize
from scipy.sparse import eye, vstack, spmatrix
from scipy.ndimage import zoom
import numpy as np
from typing import Any, List, Literal, Optional, Union
import scipy.sparse as sp
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from sklearn.preprocessing import OneHotEncoder
import logging
logging.basicConfig(format='%(asctime)s %(message)s')

COLORS = [
    [0.9019607843137255, 0.09803921568627451, 0.29411764705882354],
    [0.23529411764705882, 0.7058823529411765, 0.29411764705882354],
    [1.0, 0.8823529411764706, 0.09803921568627451],
    [0.2627450980392157, 0.38823529411764707, 0.8470588235294118],
    [0.9607843137254902, 0.5098039215686274, 0.19215686274509805],
    [0.5686274509803921, 0.11764705882352941, 0.7058823529411765],
    [0.27450980392156865, 0.9411764705882353, 0.9411764705882353],
    [0.9411764705882353, 0.19607843137254902, 0.9019607843137255],
    [0.7372549019607844, 0.9647058823529412, 0.047058823529411764],
    [0.9803921568627451, 0.7450980392156863, 0.7450980392156863],
    [0.0, 0.5019607843137255, 0.5019607843137255],
    [0.9019607843137255, 0.7450980392156863, 1.0],
    [0.6039215686274509, 0.38823529411764707, 0.1411764705882353],
    [1.0, 0.9803921568627451, 0.7843137254901961],
    [0.5019607843137255, 0.0, 0.0],
    [0.6666666666666666, 1.0, 0.7647058823529411],
    [0.5019607843137255, 0.5019607843137255, 0.0],
    [1.0, 0.8470588235294118, 0.6941176470588235],
    [0.0, 0.0, 0.4588235294117647],
    [0.5019607843137255, 0.5019607843137255, 0.5019607843137255],
    [1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0],
]

COLORS = [[int(255 * v) for v in RGB] for RGB in COLORS]




def polygons2json(polygons, cluster_class, cluster_names, colors=None):
    jsonROIs = []
    for i, polygon in enumerate(polygons):
        name = cluster_names[i]
        jsonROIs.append(
            {
                "type": "Feature",
                "geometry": {"type": "MultiPolygon", "coordinates": []},
                "properties": {
                    "name": name,
                    "classification": {"name": cluster_class},
                    "color": colors[i] if colors is not None else [255, 0, 0],
                    "isLocked": False,
                },
            }
        )
        jsonROIs[-1]["geometry"]["coordinates"].append(polygon)
    return jsonROIs

def binary_mask_to_polygon(binary_mask: np.ndarray, tolerance: float=0, offset: float=None, scale: float=None):
    """Converts a binary mask to COCO polygon representation
    Args:
        binary_mask: a 2D binary numpy array where '1's represent the object
        tolerance: Maximum distance from original points of polygon to approximated
            polygonal chain. If tolerance is 0, the original coordinate array is returned.
    """
    from skimage.measure import approximate_polygon, find_contours

    polygons = []
    # pad mask to close contours of shapes which start and end at an edge
    binary_mask = zoom(binary_mask, 3, order=0, grid_mode=True)
    padded_binary_mask = np.pad(
        binary_mask, pad_width=1, mode="constant", constant_values=0
    )
    contours = find_contours(padded_binary_mask, 0.5)
    contours = [c-1 for c in contours]
    #contours = np.subtract(contours, 1)
    for contour in contours:
        contour = approximate_polygon(contour, tolerance)
        if len(contour) < 3:
            continue
        contour = contour / 3
        contour = np.rint(contour)
        if scale is not None:
            contour = contour * scale
        if offset is not None:
            contour = contour + offset  # .ravel().tolist()

        # after padding and subtracting 1 we may get -0.5 points in our segmentation
        polygons.append(contour.tolist())

    return polygons

def labelmask2geojson(
    labelmask: np.ndarray,
    region_name:str="My regions",
    scale:float=1.0,
    offset:float = 0
):
    from skimage.measure import regionprops

    nclusters = np.max(labelmask)
    colors = [COLORS[k % len(COLORS)] for k in range(nclusters)]

    # Make JSON
    polygons = []
    cluster_names = [f"Region {l+1}" for l in np.arange(1,nclusters+1)]
    for index, region in enumerate(regionprops(labelmask)):
        # take regions with large enough areas
        contours = binary_mask_to_polygon(
            region.image,
            offset=scale*np.array(region.bbox[0:2]) + offset,
            scale=scale
        )
        polygons.append(contours)
    json = polygons2json(polygons, region_name, cluster_names, colors=colors)
    return json


def map2numeric(data: np.ndarray) -> np.ndarray:
    map2numeric = {k : i for i,k in enumerate(np.unique(data))}
    return np.array([map2numeric[v] for v in data])


def create_features(xy: np.ndarray, labels: np.ndarray, unique_labels:np.ndarray, sigma: float, bin_width: Union[float, str, None], min_genes_per_bin: int):
    if isinstance(bin_width, str):
        if bin_width == 'auto':
            bin_width = sigma/3.0

    grid_props = {}
    # Compute binning matrix
    if bin_width is not None:
        B, grid_props = spatial_binning_matrix(xy, box_width=bin_width, return_grid_props=True)
    else:
        B = eye(len(xy))
    B = B.astype('float32')

    if bin_width is not None:
        x = grid_props['grid_coords'][0] / grid_props['grid_scale'] + grid_props['grid_offset'][0]
        y = grid_props['grid_coords'][1] / grid_props['grid_scale'] + grid_props['grid_offset'][1]
        xy = np.vstack((x,y)).T


    # Create attribute matrix (ngenes x nuniques)
    attributes, _ = attribute_matrix(labels, unique_labels)
    attributes = attributes.astype('bool')

    features = kde_per_label(xy, B @ attributes, sigma)

    # Compute bin size
    bin_size = features.sum(axis=1).A.flatten()
    good_bins = bin_size >= min_genes_per_bin
    norms_r =1.0 / features.sum(axis=1)
    norms_r[np.isinf(norms_r)] = .0
    features = features.multiply(norms_r).tocsr()

    return dict(
        features=features,
        xy_bin=xy,
        norms=norms_r.A.flatten(),
        grid_props=grid_props,
        good_bins=good_bins,
        back_map=B.T.nonzero()[1]
    )

def predict(kmeans_model, features: spmatrix, good_bins: np.ndarray, back_map: np.ndarray):
    clusters = np.zeros(features.shape[0], dtype='int') - 1
    clusters[good_bins] = kmeans_model.predict(features[good_bins,:])
    return clusters[back_map], clusters

def points2regions(
        xy: np.ndarray, 
        gene_labels: np.ndarray, 
        sigma: float, 
        n_clusters: int,
        bin_width: Union[float, str, None] = 'auto', 
        min_genes_per_bin:int = 1, 
        groupids: Union[np.ndarray,None] = None, 
        convert_to_geojson: bool = False, seed:int=42, 
        region_name:str="My regions", 
        return_anndata: bool = False
    ) -> Any:
    """
    Cluster categorical points with spatial locations.

    This function takes spatial point data along with their associated labels (e.g., genes from an in situ sequencing experiment),
    performs clustering on the points, and assigns each point to a specific cluster based on the clustering results.

    Parameters:
    -----------
    xy : np.ndarray
        A 2D numpy array containing the x-y coordinates of the data points.
    gene_labels : np.ndarray
        A 1D numpy array containing the gene label associated with each data point.
    sigma : float
        The bandwidth parameter for the kernel density estimation used in clustering.
        This parameter reflects the spatial resolution of the clustering.
    n_clusters : int
        The number of clusters to generate using K-Means clustering.
    bin_width : Union[float, str, None], optional (default='auto')
        The width of the bins in which to group the points for clustering. If 'auto', the bin width is calculated as sigma / 3.0.
        If set as None, no binning is performed. This can be slow for larger sigmas.
    min_genes_per_bin : int, optional (default=1)
        The minimum number of genes required per bin for it to be considered valid. Default 1.
    groupids : Union[np.ndarray, None], optional (default=None)
        An optional array containing group IDs for each point. If provided, features are computed separately for each group.
        The whole collection of features are then clustered together. 
    convert_to_geojson : bool, optional (default=False)
        Whether to convert the clustering results to GeoJSON format, which represents the clustered regions as polygons.
    seed : int, optional (default=42)
        The random seed used for K-Means clustering.
    region_name : str, optional (default="My regions")
        The name to assign to the regions generated by clustering.
    return_anndata : bool, optional (default=True)
        Whether to return the results in the form of an AnnData object, which can store both clustering results and spatial coordinates.

    Returns:
    --------
    output : AnnData or np.ndarray
        If 'return_anndata' is True, returns an AnnData object containing clustering results and spatial coordinates. If False, returns an ndarray containing cluster assignments.
    
    Example:
    --------
    >>> import pandas as pd
    >>> data = pd.read_csv('result.csv')
    >>> xy = data[['X', 'Y']].to_numpy()
    >>> genes = data['Genes'].to_numpy()
    >>> sigma = 500
    >>> n_clusters = 8
    >>> clusters = points2regions(xy, genes, sigma, n_clusters)
    """

    
    if not isinstance(xy, np.ndarray):
        raise ValueError("Input 'xy' must be a NumPy array.")
    
    if xy.shape[1] != 2:
        raise ValueError("Input 'xy' must be a 2D numpy array with shape (n, 2), where n is the number of data points.")
    
    if not isinstance(gene_labels, np.ndarray):
        raise ValueError("Input 'gene_labels' must be a numpy array.")
    
    if gene_labels.ndim != 1:
        raise ValueError("Input 'gene_labels' must be a 1D numpy array.")
    
    if not isinstance(sigma, (int, float)):
        raise ValueError("Input 'sigma' must be a numerical value.")
    
    if not isinstance(n_clusters, int):
        raise ValueError("Input 'n_clusters' must be an integer value.")
    
    if bin_width is not None and not isinstance(bin_width, (float, str)):
        raise ValueError("Input 'bin_width' must be a float, string 'auto', or None.")
    
    if not isinstance(min_genes_per_bin, int):
        raise ValueError("Input 'min_genes_per_bin' must be an integer value.")
    
    if groupids is not None and not isinstance(groupids, np.ndarray):
        raise ValueError("Input 'groupids' must be a numpy array or None.")
    
    if not isinstance(convert_to_geojson, bool):
        raise ValueError("Input 'convert_to_geojson' must be a boolean value.")
    
    if not isinstance(seed, int):
        raise ValueError("Input 'seed' must be an integer value.")
    
    if not isinstance(region_name, str):
        raise ValueError("Input 'region_name' must be a string.")
    
    if not isinstance(return_anndata, bool):
        raise ValueError("Input 'return_anndata' must be a boolean value.")
    
    
    print (
        "xy",xy, "labels", gene_labels,"sigma",sigma, "n_clusters", n_clusters, "bin_width", bin_width, "min_genes_per_bin", min_genes_per_bin, "library_id_column", groupids, "convert_to_geojson", convert_to_geojson, "seed", seed
    )
    xy = np.array(xy, dtype="float32")

    # Iterate data by library ids
    if groupids is not None:
        unique_library_id = np.unique(groupids)
        iterdata = [
            (lib_id, (
                xy[groupids==lib_id,:],
                gene_labels[groupids==lib_id]
            )) for lib_id in unique_library_id
        ]
        get_slice = lambda library_id, data: data == library_id
    else:
        iterdata = [('id', (xy, gene_labels))]
        get_slice = lambda library_id, data: np.ones(len(data), dtype='bool')


    unique_genes = np.unique(gene_labels)
    results = {
        library_id : create_features(
            xy_slice,
            labels_slice,
            unique_genes,
            sigma,
            bin_width,
            min_genes_per_bin
        )
        for library_id, (xy_slice, labels_slice) in iterdata
    }

    # Create train features
    X_train = vstack([
        r['features'][r['good_bins']] for r in results.values()
    ])

    # Train K-Means
    logging.warning("Training KMeans")
    kmeans = KMeans(n_clusters=n_clusters, n_init='auto', random_state=seed)
    logging.warning("Fiting KMeans")
    kmeans = kmeans.fit(X_train)
    logging.warning("KMeans done")

    # Predict
    logging.warning("Predict")
    for library_id, result_dict in results.items():
        cluster_per_gene, cluster_per_bin = predict(kmeans, result_dict['features'], result_dict['good_bins'], result_dict['back_map'])
        results[library_id]['cluster_per_gene'] = cluster_per_gene
        results[library_id]['cluster_per_bin'] = cluster_per_bin

    # Add clusters to dataframe
    logging.warning("Adding clusters to dataframe")
    clusters_per_gene = np.zeros(len(xy), dtype='int')
    for library_id in results.keys():
        if groupids is not None:
            library_id_slice_ind = get_slice(library_id, groupids)
        else:
            library_id_slice_ind = get_slice(library_id, xy)
        clusters_per_gene[library_id_slice_ind] = results[library_id]['cluster_per_gene']

    if return_anndata:
        logging.warning("Exporting dataframe")
        # Create an adata object
        import anndata
        import pandas as pd

        # Get position of bins for each group (library id)
        xy_bin = np.vstack([
            r['xy_bin'][r['good_bins']] for r in results.values()
        ])

        # Get labels of bins for each group (library id)
        labels_bin = np.hstack([
            r['cluster_per_bin'][r['good_bins']] for r in results.values()
        ])

        obs = {}
        obs['points2regions'] = labels_bin
        if len(results) > 1:
            obs['groupid'] =  np.hstack([[id]*len(r['cluster_per_bin']) for id, r in results.items()])


        # Multiply back features with the norm
        norms = 1.0 / np.hstack([r['norms'][r['good_bins']] for r in results.values()])
        norms[np.isinf(norms)] = 0
        norms = norms.reshape((-1,1))

        adata = anndata.AnnData(
            X=X_train.multiply(norms).tocsc(),
            obsm={
                'spatial' : xy_bin,
            },
            obs=obs,
            var=pd.DataFrame(index=unique_genes)
        )

        adata.obs['points2regions'] = adata.obs['points2regions'].astype('category')


        if len(results) > 1:
            adata.obs['groupid'] = adata.obs['groupid'].astype('category')

        # Create reads dataframe
        reads = {}
        reads['x'] = xy[:,0]
        reads['y'] = xy[:,1]
        reads['labels'] = gene_labels
        reads['points2regions'] = clusters_per_gene

        if groupids is not None:
            reads['groupid'] = groupids

        reads = pd.DataFrame(reads).reset_index(drop=True)
        reads['labels'] = reads['labels'].astype('category')
        reads['points2regions'] = reads['points2regions'].astype('category')
        if groupids is not None:
            reads['groupid'] = reads['groupid'].astype('category')
        adata.uns['reads'] = reads
        output = adata

    else:
        output = clusters_per_gene      

    if convert_to_geojson:
        geojsons = []
        for result in results.values():
            grid_props = result['grid_props']
            clusters = result['cluster_per_bin']
            label_mask = np.zeros(grid_props['grid_size'], dtype='uint8')
            label_mask[tuple(ind for ind in grid_props['grid_coords'])] = clusters
            label_mask = label_mask
            geojson = labelmask2geojson(label_mask, region_name=region_name, scale=1.0/grid_props['grid_scale'], offset=grid_props['grid_offset'])
            geojsons.append(geojson)
        return (output, geojsons)
    else:
        return output


def connectivity_matrix(
    xy: np.ndarray,
    method="knn",
    k: int = 5,
    r: Optional[float] = None,
    include_self: bool = False,
) -> sp.spmatrix:
    """
    Compute the connectivity matrix of a dataset based on either k-NN or radius search.

    Parameters
    ----------
    xy : np.ndarray
        The input dataset, where each row is a sample point.
    method : str, optional (default='knn')
        The method to use for computing the connectivity.
        Can be either 'knn' for k-nearest-neighbors or 'radius' for radius search.
    k : int, optional (default=5)
        The number of nearest neighbors to use when method='knn'.
    r : float, optional (default=None)
        The radius to use when method='radius'.
    include_self : bool, optional (default=False)
        If the matrix should contain self connectivities.

    Returns
    -------
    A : sp.spmatrix
        The connectivity matrix, with ones in the positions where two points are
            connected.
    """
    if method == "knn":
        A = kneighbors_graph(xy, k, include_self=include_self).astype('bool')
    else:
        A = radius_neighbors_graph(xy, r, include_self=include_self).astype('bool')
    return A


def attribute_matrix(
    cat: np.ndarray,
    unique_cat: Union[np.ndarray, Literal["auto"]] = "auto",
    return_encoder: bool = False,
):
    """
    Compute the attribute matrix from categorical data, based on one-hot encoding.

    Parameters
    ----------
    cat : np.ndarray
        The categorical data, where each row is a sample and each column is a feature.
    unique_cat : np.ndarray
        Unique categorical data used to setup up the encoder. If "auto", unique
        categories are automatically determined from cat.
    return_encoder : bool, optional (default=False)
        Whether to return the encoder object, in addition to the attribute matrix and
        categories list.

    Returns
    -------
    y : sp.spmatrix
        The attribute matrix, in sparse one-hot encoding format.
    categories : list
        The categories present in the data, as determined by the encoder.
    encoder : OneHotEncoder
        The encoder object, only returned if \`return_encoder\` is True.
    """
    X = np.array(cat).reshape((-1, 1))
    if not isinstance(unique_cat, str):
        unique_cat_list: Union[List[np.ndarray], Literal["auto"]] = [
            np.array(unique_cat)
        ]
    elif unique_cat == "auto":
        unique_cat_list = "auto"
    else:
        raise ValueError("\`unique_cat\` must be a numpy array or the string \`auto\`.")
    encoder = OneHotEncoder(
        categories=unique_cat_list, sparse_output=True, handle_unknown="ignore"
    )
    encoder.fit(X)
    y = encoder.transform(X)
    categories = list(encoder.categories_[0])
    if return_encoder:
        return y, categories, encoder
    return y, categories




def spatial_binning_matrix(
    xy: np.ndarray, box_width: float, return_grid_props: bool = False
) -> sp.spmatrix:
    """
    Compute a sparse matrix that indicates which points in a point cloud fall in which
    hyper-rectangular bins.

    Parameters:
    points (numpy.ndarray): An array of shape (N, D) containing the D-dimensional
        coordinates of N points in the point cloud.
    box_width (float): The width of the bins in which to group the points.

    Returns:
    sp.spmatrix: A sparse matrix of shape (num_bins, N) where num_bins is the number of
        bins. The matrix is such that the entry (i,j) is 1 if the j-th point falls in
        the i-th bin, and 0 otherwise.

    Example:
    >>> points = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [1.5, 1.5, 1.5], [2, 2, 2]])
    >>> bin_matrix = spatial_binning_matrix(points, 1)
    >>> print(bin_matrix.toarray())
    [[1 1 0 0]
        [1 1 0 0]
        [0 0 1 1]]
    """

    # Compute shifted coordinates
    mi, ma = xy.min(axis=0, keepdims=True), xy.max(axis=0, keepdims=True)
    xys = xy - mi

    # Compute grid size
    grid = ma - mi
    grid = grid.flatten()

    # Compute bin index
    bin_ids = xys // box_width
    bin_ids = bin_ids.astype("int")
    bin_ids = tuple(x for x in bin_ids.T)

    # Compute grid size in indices
    size = grid // box_width + 1
    size = tuple(x for x in size.astype("int"))

    # Convert bin_ids to integers
    linear_ind = np.ravel_multi_index(bin_ids, size)

    # Create a matrix indicating which markers fall in what bin
    bin_matrix, linear_unique_bin_ids = attribute_matrix(linear_ind)


    bin_matrix = bin_matrix.T

    if return_grid_props:
        sub_unique_bin_ids = np.unravel_index(linear_unique_bin_ids, size)
        grid_props = dict(
            grid_coords=sub_unique_bin_ids,
            grid_size=size,
            grid_offset=mi.flatten(),
            grid_scale=1.0/box_width
        )

    return (bin_matrix, grid_props) if return_grid_props else bin_matrix

def kde_per_label(xy: np.ndarray, features: sp.spmatrix, sigma: float, return_neighbors: bool = False):
    """
    Computes the kernel density estimation (KDE) for each label in \`labels\`, using the
    data points in \`xy\` as inputs. Returns the KDE values as an attribute matrix, and
    the unique labels found in \`labels\`.

    Parameters:
    -----------
    xy : numpy.ndarray
        A 2D numpy array of shape (n, 2) containing the x-y coordinates of the data
        points.
    features : sp.spmatrix
        Features that are to be blured using KDE
    sigma : float
        The standard deviation of the Gaussian kernel to use in the KDE.

    Returns:
    --------
    Tuple of two numpy.ndarray:
        - \`att\`: A 2D numpy array of shape (n_labels, n_features), where n_labels is the
            number of unique labels in \`labels\`
                  and n_features is the number of attributes (columns) in \`labels\`. Each
                  row represents the KDE values for a single label.
        - \`unique_labels\`: A 1D numpy array containing the unique labels found in
            \`labels\`.
    """
    logging.warning ("Compute connectivity matrix")
    adj = connectivity_matrix(xy, method="radius", r=2.0 * sigma, include_self=True)
    logging.warning ("Compute connectivity matrix done")
    row, col = adj.nonzero()
    d2 = (xy[row,0] - xy[col,0])**2
    d2 = d2 + (xy[row,1] - xy[col,1])**2
    d2 = np.exp(-d2 / (2 * sigma * sigma))
    logging.warning ("Compute csr_matrix")
    aff = sp.csr_matrix((d2, (row, col)), shape=adj.shape, dtype='float32')
    logging.warning ("Compute csr_matrix done")
    if not return_neighbors:
        return aff @ features
    else:
        return aff @ features, adj




if __name__ == '__main__':


    import pandas as pd
    data = pd.read_csv('result.csv')
    #data1['ID'] = 'asd'
    #data2 = pd.read_csv('result.csv').query("X > 5000 and X < 10000 and Y > 5000 and Y <= 10000")
    #data2['ID'] = 'basd'
    #data = pd.concat((data1, data2), axis=0)
    xy = data[['X', 'Y']].to_numpy()
    genes = data['Genes'].to_numpy()
    sigma = 10
    bin_width = 50.0
    min_genes_per_bin = 10
    n_clusters = 8
    adata = points2regions(xy, genes, sigma, n_clusters, bin_width, min_genes_per_bin, return_anndata=True)
    adata.write_h5ad('test.h5ad')
    
    clusters = adata.obs.copy()
    clusters['x'] = adata.obsm['spatial'][:,0]
    clusters['y'] = adata.obsm['spatial'][:,1]
    clusters.to_csv('clusters.csv', index=False)
    adata.uns['reads'].to_csv('reads.csv', index=False)
