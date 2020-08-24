import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/paul15.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures
adata = sc.datasets.paul15()
adata.X = adata.X.astype('float64')  # this is not required and results will be comparable without it
sc.pp.recipe_zheng17(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)

pl.figure(figsize=(12, 12))
sc.pl.draw_graph(adata, color='paul15_clusters') #legend_loc='on data'
pl.savefig("paul15_clusters")



#paga
sc.tl.louvain(adata, resolution=1.0)
sc.tl.paga(adata, groups='louvain')

pl.figure(figsize=(12, 12))
sc.pl.paga(adata, color=['louvain', 'Hba-a2', 'Elane', 'Irf8'])
pl.savefig("paul15_paga.pdf")



