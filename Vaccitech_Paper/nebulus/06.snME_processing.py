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
    binarize_matrix
from ALLCools.plot import *
from ALLCools.mcds import MCDS

mcds_path = '/Genomics/ogtr04/wat2/snME/data/vaccitech_flu_snME_final.mcds/'

# PC cutoff
pc_cutoff = 0.1

resolution = 1

mcds = MCDS.open(mcds_path, var_dim='chrom5k')

mcad = mcds.get_score_adata(mc_type='CGN', quant_type='hypo-score')

binarize_matrix(mcad, cutoff=0.95)

filter_regions(mcad)

mcad

# by default we save the results in adata.obsm['X_pca'] which is the scanpy defaults in many following functions
# But this matrix is not calculated by PCA
lsi(mcad, algorithm='arpack', obsm='X_pca')

# choose significant components
n_components = significant_pc_test(mcad, p_cutoff=pc_cutoff, update=True)

sc.pp.neighbors(mcad)
sc.tl.leiden(mcad, resolution=resolution)

try:
    sc.tl.paga(mcad, groups='leiden')
    sc.pl.paga(mcad, plot=False)
    sc.tl.umap(mcad, init_pos='paga')
except:
    sc.tl.umap(mcad)
    
fig, ax = plt.subplots(figsize=(4, 4), dpi=250)
_ = categorical_scatter(data=mcad, ax=ax, coord_base='umap', hue='leiden', show_legend=True)

interactive_scatter(data=mcad, hue='leiden', coord_base='umap')

mcad.write_h5ad(f'adata.with_coords.mcad')

mcad