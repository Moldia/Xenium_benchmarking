""" SpaGE [1]
@author: Tamim Abdelaal
This function integrates two single-cell datasets, spatial and scRNA-seq, and 
enhance the spatial data by predicting the expression of the spatially 
unmeasured genes from the scRNA-seq data.
The integration is performed using the domain adaption method PRECISE [2]
	
References
----------
    [1] Abdelaal T., Mourragui S., Mahfouz A., Reiders M.J.T. (2020)
    SpaGE: Spatial Gene Enhancement using scRNA-seq
    [2] Mourragui S., Loog M., Reinders M.J.T., Wessels L.F.A. (2019)
    PRECISE: A domain adaptation approach to transfer predictors of drug response
    from pre-clinical models to tumors
"""

import numpy as np
import pandas as pd
import scipy.stats as st
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

def SpaGE(Spatial_data,RNA_data,n_pv,genes_to_predict=None):
    """
        @author: Tamim Abdelaal
        This function integrates two single-cell datasets, spatial and scRNA-seq, 
        and enhance the spatial data by predicting the expression of the spatially 
        unmeasured genes from the scRNA-seq data.
        
        Args:  
            Spatial_data : Dataframe
                Normalized Spatial data matrix (cells X genes).

            RNA_data : Dataframe
                Normalized scRNA-seq data matrix (cells X genes).

            n_pv : int
                Number of principal vectors to find from the independently computed
                principal components, and used to align both datasets. This should
                be <= number of shared genes between the two datasets.

            genes_to_predict : str array 
                list of gene names missing from the spatial data, to be predicted 
                from the scRNA-seq data. Default is the set of different genes 
                (columns) between scRNA-seq and spatial data.
            
        results:  
            Imp_Genes: Dataframe
                Matrix containing the predicted gene expressions for the spatial 
                cells. Rows are equal to the number of spatial data rows (cells), 
                and columns are equal to genes_to_predict,  .
    """
    
    if genes_to_predict is SpaGE.__defaults__[0]:
        genes_to_predict = np.setdiff1d(RNA_data.columns,Spatial_data.columns)
        
    RNA_data_scaled = pd.DataFrame(data=st.zscore(RNA_data,axis=0),
                                   index = RNA_data.index,columns=RNA_data.columns)
    Spatial_data_scaled = pd.DataFrame(data=st.zscore(Spatial_data,axis=0),
                                   index = Spatial_data.index,columns=Spatial_data.columns)
    Common_data = RNA_data_scaled[np.intersect1d(Spatial_data_scaled.columns,RNA_data_scaled.columns)]

    Y_train = RNA_data[genes_to_predict]
    
    Imp_Genes = pd.DataFrame(np.zeros((Spatial_data.shape[0],len(genes_to_predict))),
                                 columns=genes_to_predict)
    
    pv_Spatial_RNA = PVComputation(
            n_factors = n_pv,
            n_pv = n_pv,
            dim_reduction = 'pca',
            dim_reduction_target = 'pca'
    )
    
    pv_Spatial_RNA.fit(Common_data,Spatial_data_scaled[Common_data.columns])
    
    S = pv_Spatial_RNA.source_components_.T
        
    Effective_n_pv = sum(np.diag(pv_Spatial_RNA.cosine_similarity_matrix_) > 0.3)
    S = S[:,0:Effective_n_pv]
    
    Common_data_projected = Common_data.dot(S)
    Spatial_data_projected = Spatial_data_scaled[Common_data.columns].dot(S)
        
    nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',
                            metric = 'cosine').fit(Common_data_projected)
    distances, indices = nbrs.kneighbors(Spatial_data_projected)
    
    for j in range(0,Spatial_data.shape[0]):
    
        weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
        weights = weights/(len(weights)-1)
        Imp_Genes.iloc[j,:] = np.dot(weights,Y_train.iloc[indices[j,:][distances[j,:] < 1]])
        
    return Imp_Genes

""" Dimensionality Reduction
@author: Soufiane Mourragui
This module extracts the domain-specific factors from the high-dimensional omics
dataset. Several methods are here implemented and they can be directly
called from string name in main method method. All the methods
use scikit-learn implementation.
Notes
-------
	-
	
References
----------
	[1] Pedregosa, Fabian, et al. (2011) Scikit-learn: Machine learning in Python.
	Journal of Machine Learning Research
"""

import numpy as np
from sklearn.decomposition import PCA, FastICA, FactorAnalysis, NMF, SparsePCA
from sklearn.cross_decomposition import PLSRegression


def process_dim_reduction(method='pca', n_dim=10):
    """ Default linear dimensionality reduction method. For each method, return a
    BaseEstimator instance corresponding to the method given as input.

    Args:
        method: str, default to 'pca'
            Method used for dimensionality reduction.
            Implemented: 'pca', 'ica', 'fa' (Factor Analysis), 
            'nmf' (Non-negative matrix factorisation), 'sparsepca' (Sparse PCA).
        
        n_dim: int, default to 10
            Number of domain-specific factors to compute.

    results:
        Classifier, i.e. BaseEstimator instance
    """

    if method.lower() == 'pca':
        clf = PCA(n_components=n_dim)

    elif method.lower() == 'ica':
        print('ICA')
        clf = FastICA(n_components=n_dim)

    elif method.lower() == 'fa':
        clf = FactorAnalysis(n_components=n_dim)

    elif method.lower() == 'nmf':
        clf = NMF(n_components=n_dim)

    elif method.lower() == 'sparsepca':
        clf = SparsePCA(n_components=n_dim, alpha=10., tol=1e-4, verbose=10, n_jobs=1)

    elif method.lower() == 'pls':
        clf = PLS(n_components=n_dim)
		
    else:
        raise NameError('%s is not an implemented method'%(method))

    return clf


class PLS():
    """
    Implement PLS to make it compliant with the other dimensionality
    reduction methodology.
    (Simple class rewritting).
    """
    def __init__(self, n_components=10):
        self.clf = PLSRegression(n_components)

    def get_components_(self):
        return self.clf.x_weights_.transpose()

    def set_components_(self, x):
        pass

    components_ = property(get_components_, set_components_)

    def fit(self, X, y):
        self.clf.fit(X,y)
        return self

    def transform(self, X):
        return self.clf.transform(X)

    def predict(self, X):
        return self.clf.predict(X)
    
    
    
""" Principal Vectors
@author: Soufiane Mourragui
This module computes the principal vectors from two datasets, i.e.:
- perform linear dimensionality reduction independently for both dataset, resulting
in set of domain-specific factors.
- find the common factors using principal vectors [1]
This result in set of pairs of vectors. Each pair has one vector from the source and one
from the target. For each pair, a similarity score (cosine similarity) can be computed
between the principal vectors and the pairs are naturally ordered by decreasing order
of this similarity measure.
Example
-------
    Examples are given in the vignettes.
Notes
-------
	Examples are given in the vignette
	
References
----------
	[1] Golub, G.H. and Van Loan, C.F., 2012. "Matrix computations" (Vol. 3). JHU Press.
	[2] Mourragui, S., Loog, M., Reinders, M.J.T., Wessels, L.F.A. (2019)
    PRECISE: A domain adaptation approach to transfer predictors of drug response
    from pre-clinical models to tumors
"""

import numpy as np
import pandas as pd
import scipy
from pathlib import Path
from sklearn.preprocessing import normalize


class PVComputation:
    """
    Attributes:
        n_factors: int
            Number of domain-specific factors to compute.

        n_pv: int
            Number of principal vectors.

        dim_reduction_method_source: str
            Dimensionality reduction method used for source data.

        dim_reduction_target: str
            Dimensionality reduction method used for source data.

        source_components_ : numpy.ndarray, shape (n_pv, n_features)
            Loadings of the source principal vectors ranked by similarity to the
            target. Components are in the row.

        source_explained_variance_ratio_: numpy.ndarray, shape (n_pv)
            Explained variance of the source on each source principal vector.

        target_components_ : numpy.ndarray, shape (n_pv, n_features)
            Loadings of the target principal vectors ranked by similarity to the
            source. Components are in the row.

        target_explained_variance_ratio_: numpy.ndarray, shape (n_pv)
            Explained variance of the target on each target principal vector.

        cosine_similarity_matrix_: numpy.ndarray, shape (n_pv, n_pv)
            Scalar product between the source and the target principal vectors. Source
            principal vectors are in the rows while target's are in the columns. If
            the domain adaptation is sensible, a diagonal matrix should be obtained.
    """

    def __init__(self, n_factors,n_pv,
                dim_reduction='pca',
                dim_reduction_target=None,
                project_on=0):
        """
        Args:
            n_factors : int
                Number of domain-specific factors to extract from the data (e.g. using PCA, ICA).

            n_pv : int
                Number of principal vectors to find from the independently computed factors.

            dim_reduction : str, default to 'pca' 
                Dimensionality reduction method for the source data,
                i.e. 'pca', 'ica', 'nmf', 'fa', 'sparsepca', pls'.

            dim_reduction_target : str, default to None 
                Dimensionality reduction method for the target data,
                i.e. 'pca', 'ica', 'nmf', 'fa', 'sparsepca', pls'. If None, set to dim_reduction.

            project_on: int or bool, default to 0
                Where data should be projected on. 0 means source PVs, -1 means target PVs and 1 means
                both PVs.
        """

        self.n_factors = n_factors
        self.n_pv = n_pv
        self.dim_reduction_method_source = dim_reduction
        self.dim_reduction_method_target = dim_reduction_target or dim_reduction
        self.dim_reduction_source = self._process_dim_reduction(self.dim_reduction_method_source)
        self.dim_reduction_target = self._process_dim_reduction(self.dim_reduction_method_target)

        self.source_components_ = None
        self.source_explained_variance_ratio_ = None
        self.target_components_ = None
        self.target_explained_variance_ratio_ = None
        self.cosine_similarity_matrix_ = None

    def _process_dim_reduction(self, dim_reduction):
        if type(dim_reduction) == str:
            return process_dim_reduction(method=dim_reduction, n_dim=self.n_factors)
        else:
            return dim_reduction

    def fit(self, X_source, X_target, y_source=None):
        """ Compute the common factors between two set of data.
        IMPORTANT: Same genes have to be given for source and target, and in same order

        Args:
            X_source : np.ndarray, shape (n_components, n_genes)
                Source dataset.

            X_target : np.ndarray, shape (n_components, n_genes)
                Target dataset.

            y_source : np.ndarray, shape (n_components, 1) (optional, default to None)
                Eventual output, in case one wants to give ouput (for instance PLS).

        results:
            self: returns an instance of self.
        """

        # Compute factors independently for source and target. Orthogonalize the basis
        Ps = self.dim_reduction_source.fit(X_source, y_source).components_
        Ps = scipy.linalg.orth(Ps.transpose()).transpose()

        Pt = self.dim_reduction_target.fit(X_target, y_source).components_
        Pt = scipy.linalg.orth(Pt.transpose()).transpose()

        # Compute the principal factors
        self.compute_principal_vectors(Ps, Pt)

        # Compute variance explained
        self.source_explained_variance_ratio_ = np.var(self.source_components_.dot(X_source.transpose()), axis=1)/\
                                                np.sum(np.var(X_source), axis=0)
        self.target_explained_variance_ratio_ = np.var(self.target_components_.dot(X_target.transpose()), axis=1)/\
                                                np.sum(np.var(X_target), axis=0)

        return self

    def compute_principal_vectors(self, source_factors, target_factors):
        """Compute the principal vectors between the already computed set of domain-specific
        factors, using approach presented in [1,2].
        IMPORTANT: Same genes have to be given for source and target, and in same order

        Args:
            source_factors: np.ndarray, shape (n_components, n_genes)
                Source domain-specific factors.

            target_factors: np.ndarray, shape (n_components, n_genes)
                Target domain-specific factors.

        results:
            self: returns an instance of self.
        """

        # Find principal vectors using SVD
        u,sigma,v = np.linalg.svd(source_factors.dot(target_factors.transpose()))
        self.source_components_ = u.transpose().dot(source_factors)[:self.n_pv]
        self.target_components_ = v.dot(target_factors)[:self.n_pv]
        # Normalize to make sure that vectors are unitary
        self.source_components_ = normalize(self.source_components_, axis = 1)
        self.target_components_ = normalize(self.target_components_, axis = 1)

        # Compute cosine similarity matrix
        self.initial_cosine_similarity_matrix_ = source_factors.dot(target_factors.transpose())
        self.cosine_similarity_matrix_ = self.source_components_.dot(self.target_components_.transpose())

        # Compute angles
        self.angles_ = np.arccos(np.diag(self.cosine_similarity_matrix_))

        return self


    def transform(self, X, project_on=None):
        """ Projects data onto principal vectors.
        
        Args:
            X : numpy.ndarray, shape (n_samples, n_genes)
                Data to project.

            project_on: int or bool, default to None
                Where data should be projected on. 0 means source PVs, -1 means target PVs and 1 means
                both PVs. If None, set to class instance value.

        results:
            Projected data as a numpy.ndarray of shape (n_samples, n_factors).
        """

        project_on = project_on or self.project_on

        # Project on source
        if project_on == 0:
            return X.dot(self.source_components_.transpose())

        # Project on target
        elif project_on == -1:
            return X.dot(self.target_components_.transpose())

        # Project on both
        elif project_on == 1:
            return X.dot(np.concatenate([self.source_components_.transpose(), self.target_components_.transpose()]))

        else:
            raise ValueError('project_on should be 0 (source), -1 (target) or 1 (both). %s not correct value'%(project_on))



def leave_one_out_validation(adata,sc_adata,genes:list):
    '''Function to validate the imputation of genes using SpaGe'''
    RNA_data=sc_adata.to_df()
    RNA_meta=sc_adata.obs
    spatial_data=adata.to_df()
    spatial_meta=adata.obs
    
    Correlations = pd.Series(index = genes)
    plt.style.use('dark_background')
    spatialcoords=pd.DataFrame(adata.obsm['spatial'])
    allimp=[]
    for i in genes:
        Imp_Genes = SpaGE(spatial_data.drop(i,axis=1),RNA_data,n_pv=30,
                               genes_to_predict = [i])
        Correlations[i] = st.spearmanr(spatial_data[i],Imp_Genes[i])[0]

        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.axis('off')
        cmap = spatial_data[i]
        cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
        ax1.scatter(spatialcoords.iloc[:,0],spatialcoords.iloc[:,1],s=1,c=cmap)
        ax1.set_title('Measured ' + i, fontsize = 12)
        ax1.set_ylabel(i)
        ax2.axis('off')
        cmap = Imp_Genes[i]
        cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
        ax2.scatter(spatialcoords.iloc[:,0],spatialcoords.iloc[:,1],s=1,c=cmap)
        ax2.set_title('Predicted ' + i, fontsize = 12)
        allimp.append(Imp_Genes)
    return pd.concat(allimp)

def gene_imputation(adata,sc_adata,new_genes:list):
    '''Function to impute genes using SpaGe'''
    RNA_data=sc_adata.to_df()
    RNA_meta=sc_adata.obs
    spatial_data=adata.to_df()
    spatial_meta=adata.obs
    Correlations = pd.Series(index = new_genes)
    plt.style.use('dark_background')
    spatialcoords=pd.DataFrame(adata.obsm['spatial'])
    Imp_New_Genes = SpaGE(spatial_data,RNA_data,n_pv=30,genes_to_predict = new_genes)
    return Imp_New_Genes