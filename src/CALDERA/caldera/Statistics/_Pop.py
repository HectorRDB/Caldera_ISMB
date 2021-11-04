# Data structure:
# n samples and N nodes.
# Input:
#  - A matrix with the nodes by sample pattern of presence absence. It's an
#    N by n binary matrix.
# Output:
# - A txt file called pop.txt with one column. Each cell is a categorical
#   variable containing sample assignation of the associated sample.

# TODO
# - Need some parallelization of the k-means algorithm (need to see how to pass
#   that as on option) # njobs = -1 ?
# - Check how to find proper value of k (silhouette width? Or maybe depending on
#   phenotype structure?). # I would do silhouette also
# - Test stability w.r.t. k and D
# - Convert to a script with arguments so we can change the parameters.

import sklearn.preprocessing as prep
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import numpy as np

def assignPops(U, k, cutoff, seed = 7361):
    """
    Given a sample by node matrix, perform k-means on the matrix projected on a reduced space through PCA.
    
    Args
        U: The sample by nodes numpy array
        k: number of clusters used for k-means
        cutoff: determine how many PCs to use. We keep the top D PCs that explain
        this proportion of the total variance.
        seed: for reproducibiliy purposes
    Return
        Pop assignment using PCA and Kmeans
    """
    # Compute the decomposition
    scaledU = prep.scale(U, axis=0, with_mean=True, with_std=True)
    pca = PCA(svd_solver='full')
    pca.fit(scaledU)
    D = np.where(np.cumsum(pca.explained_variance_ratio_) > cutoff)[0]
    reducedDim = pca.transform(scaledU)[:,0:max(D)]  
    
    # Run the k-mean algorithm
    kmeans = KMeans(n_clusters = k, random_state = seed)
    pop = kmeans.fit_predict(reducedDim)

    return pop
