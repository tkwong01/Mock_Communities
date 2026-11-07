# /Users/tkwong01/Desktop/MockGMM/timeseries_counts_custom_2sp_nonzero_ORIGINAL.tsv
# First column is Reference_ID, Length, T1, T2, T3 up to T30
import numpy as np
from sklearn.decomposition import MiniBatchDictionaryLearning

def forward_orthogonal_deviations(W):
    """
    Computes Forward Orthogonal Deviations across time (axis=1).
    W shape: (P, T)
    Returns: W_fod shape (P, T-1)
    """
    P, T = W.shape
    W_fod = np.zeros((P, T - 1))
    
    for t in range(T - 1):
        r = T - (t + 1)  # Number of future periods
        c_t = np.sqrt(r / (r + 1))
        future_mean = np.mean(W[:, t + 1:], axis=1)
        W_fod[:, t] = c_t * (W[:, t] - future_mean)
        
    return W_fod

def run_metagenomic_pipeline(X_counts, M_estimated, transform_method='fod'):
    """
    X_counts: (P, T) matrix of raw gene counts (P genes, T time points)
    M_estimated: Number of hidden bacterial taxa (M) inferred from rank/BIC
    transform_method: 'fod' for Forward Orthogonal Deviations, or 'diff' for First Differences
    """
    P, T = X_counts.shape
    
    # Step 1 & 2: Pseudo-count & CLR Transformation
    eps = 1e-6
    X_pseudo = X_counts + eps
    geo_means = np.exp(np.mean(np.log(X_pseudo), axis=0))
    W = np.log(X_pseudo / geo_means) # (P, T)
    
    # Step 3: Conditional Transformation for Stationarity
    if transform_method.lower() == 'fod':
        Delta_W = forward_orthogonal_deviations(W) # (P, T - 1)
    elif transform_method.lower() == 'diff':
        Delta_W = np.diff(W, axis=1) # (P, T - 1)
    else:
        raise ValueError("transform_method must be either 'fod' or 'diff'")
    
    # Step 4: Covariance & Spectral Orthogonalization
    Sigma = np.cov(Delta_W) # (P, P)
    eigenvalues, V = np.linalg.eigh(Sigma)
    
    # Drop near-zero/negative numerical noise eigenvalues
    valid_idx = eigenvalues > 1e-8
    V_filtered = V[:, valid_idx]
    
    # Project into orthogonal coordinate system Z
    Z = np.dot(V_filtered.T, Delta_W)
    
    # Step 5: Sparse Matrix Factorization (Dictionary Learning) to find C matrix
    transformer = MiniBatchDictionaryLearning(
        n_components=M_estimated,
        alpha=0.5,
        fit_algorithm='lasso_cd',
        transform_algorithm='lasso_lars',
        random_state=42
    )
    
    S = transformer.fit_transform(Z.T).T  # (M, T-1) hidden taxon trajectories
    C = transformer.components_.T       # (K, M) weight matrix linking clusters to taxa
    
    # Step 6: Classification based on row-wise non-zero counts in C
    tolerance = 1e-4
    non_zero_counts = np.sum(np.abs(C) > tolerance, axis=1)
    
    unique_clusters = np.where(non_zero_counts == 1)[0]
    shared_clusters = np.where(non_zero_counts > 1)[0]
    background_clusters = np.where(non_zero_counts == 0)[0]
    
    return {
        "C_matrix": C,
        "hidden_taxa_trajectories": S,
        "unique_cluster_indices": unique_clusters,
        "shared_cluster_indices": shared_clusters,
        "background_cluster_indices": background_clusters
    }
