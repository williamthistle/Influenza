import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
import scanpy as sc

from ALLCools.clustering import \
    tsne, \
    significant_pc_test, \
    filter_regions, \
    remove_black_list_region, \
    lsi, \
    binarize_matrix, \
    ConsensusClustering, \
    Dendrogram, \
    get_pc_centers
from ALLCools.plot import *
from ALLCools.mcds import MCDS

# Location of MCDS file
mcds_path = '/Genomics/ogtr04/wat2/snME/data/vaccitech_flu_snME_final.mcds/'

# Parameter used for choosing PC cutoff
pc_cutoff = 0.1

# Resolution for UMAP
resolution = 1

# Open MCDS file
mcds = MCDS.open(mcds_path, var_dim='chrom5k')

# Get specific type of score we want
mcad = mcds.get_score_adata(mc_type='CGN', quant_type='hypo-score')

# Preprocessing
binarize_matrix(mcad, cutoff=0.95)

filter_regions(mcad)

# Run LSI
# by default we save the results in adata.obsm['X_pca'] which is the scanpy defaults in many following functions
# But this matrix is not calculated by PCA
lsi(mcad, algorithm='arpack', obsm='X_pca')

mcad.write_h5ad(f'adata.with_lsi.mcad')

# Subset to significant components
n_components = significant_pc_test(mcad, p_cutoff=pc_cutoff, update=True)

# Neighbors / Leiden clustering / UMAP
sc.pp.neighbors(mcad)
sc.tl.leiden(mcad, resolution=resolution)

try:
    sc.tl.paga(mcad, groups='leiden')
    sc.pl.paga(mcad, plot=False)
    sc.tl.umap(mcad, init_pos='paga')
except:
    sc.tl.umap(mcad)
    
# Plot by cluster
fig, ax = plt.subplots(figsize=(8, 6), dpi=250)
_ = categorical_scatter(data=mcad, ax=ax, coord_base='umap', hue='leiden', show_legend=True)
plt.tight_layout()
fig.savefig("snME_with_leiden_clustering.png")

mcad.write_h5ad(f'adata.with_coords.mcad')

# Add cell type annotations and write to file
metadata = pd.read_csv("snME_metadata.csv", index_col=0)
mcad.obs['CellTypeAnno'] = metadata['CellTypeAnno']
mcad.obs.to_csv('snME_clustering_results.csv.gz')

# Plot by cell type annotation
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
_ = categorical_scatter(data=mcad,
    ax=ax,
    coord_base='umap',
    hue='CellTypeAnno',
    text_anno='CellTypeAnno',
    palette='tab20',
    show_legend=True)
plt.tight_layout()
fig.savefig("snME_with_cell_types.png")


# Improve clustering via consensus clustering
# Name where cluster labels will be saved
clustering_name = 'L1'

# Basic info and read in input data
cell_meta_path = 'snME_metadata.csv'
adata_path = 'adata.with_coords.mcad'
coord_base = 'umap'

adata = anndata.read_h5ad(adata_path)

# Consensus clustering parameters
n_neighbors = 25
metric = 'euclidean'
min_cluster_size = 10
consensus_rate = 0.5
leiden_repeats = 500
leiden_resolution = 0.5
random_state = 0
n_jobs = 24
train_frac = 0.5
train_max_n = 500
max_iter = 20

# Create consensus clustering model
cc = ConsensusClustering(model=None,
                         n_neighbors=n_neighbors,
                         metric=metric,
                         min_cluster_size=min_cluster_size,
                         leiden_repeats=leiden_repeats,
                         leiden_resolution=leiden_resolution,
                         consensus_rate=consensus_rate,
                         random_state=random_state,
                         train_frac=train_frac,
                         train_max_n=train_max_n,
                         max_iter=max_iter,
                         n_jobs=n_jobs)
                         
# Run consensus clustering
cc.fit_predict(adata.obsm['X_pca'])

# Add consensus clustering label to object
adata.obs[clustering_name] = cc.label

# Plot by consensus cluster
fig, ax = plt.subplots(figsize=(8, 6), dpi=250)
_ = categorical_scatter(data=adata,
                        ax=ax,
                        hue=clustering_name,
                        coord_base=coord_base,
                        palette='tab20',
                        text_anno=clustering_name,
                        show_legend=True)
plt.tight_layout()
fig.savefig("snME_with_consensus_clustering.png")