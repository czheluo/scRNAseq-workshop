import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scprep

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')
results_file = 'write/pbmc3k.h5ad'
adata=sc.read("pbmc3k.h5ad")
id=pd.DataFrame(adata.X)
id.to_csv('pbmc3k.csv')
var=pd.DataFrame(adata.var)
var.to_csv('pbmc3k.gene.csv')
obs=pd.DataFrame(adata.obs)
obs.to_csv('pbmc3k.cell.csv')

clusters=pd.read_csv("clust_retinal_bipolar.txt",sep="\t")
adata = sc.read_10x_mtx(
    'filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)    
adata.var_names_make_unique()
adata
plt.figure(figsize=(12, 12))
sc.pl.highest_expr_genes(adata, n_top=20, )
plt.savefig("top20.png")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

plt.figure(figsize=(12, 12))
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("violin.png")
plt.savefig("violin.pdf")

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

plt.figure(figsize=(12, 12))
sc.pl.highly_variable_genes(adata)
plt.savefig("high.png")
plt.savefig("high.pdf")
adata.raw = adata

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')


sc.pl.pca(adata, color='CST3')

plt.figure(figsize=(12, 12))
sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig("pca.png")

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#tl.paga(adata)
#pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#tl.umap(adata, init_pos='paga')

sc.tl.umap(adata)
plt.figure(figsize=(12, 12))

sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])

plt.savefig("umap.png")

#Clustering the neighborhood graph
sc.tl.leiden(adata)
plt.figure(figsize=(12, 12))
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
plt.savefig("leiden.png")

adata.write(results_file)

#Finding marker genes

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

plt.figure(figsize=(12, 12))
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("markergene.png")

sc.settings.verbosity = 2  # reduce the verbosity

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

plt.figure(figsize=(12, 12))
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("markergene.wilcoxn.png")
adata.write(results_file)


sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
adata = sc.read(results_file)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

adata = sc.read(results_file)
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden')

new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
adata.rename_categories('leiden', new_cluster_names)

sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')
sc.pl.dotplot(adata, marker_genes, groupby='leiden',save='.pdf')
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90,save='.pdf')


import louvain
import scipy
from scipy import stats
import sklearn
import sklearn.cluster
import sklearn.datasets
import phate
import umap
import scprep
import graphtools as gt
import phenograph
from sklearn import cluster
from sklearn.cluster import AgglomerativeClustering

##pca

data_pca = scprep.reduce.pca(adata.X, n_components=40)

phate_op = phate.PHATE(knn=5, n_jobs=-2)
data_phate = phate_op.fit_transform(data_pca)

tic = time.time()

kmeans_clusters = sklearn.cluster.KMeans(n_clusters=8).fit_predict(data_pca)
print('Finished KMeans clustering in {:.2f} seconds'.format(time.time() - tic))

cluster_cmap = plt.cm.tab20(np.linspace(0, 1, 20))

ward = cluster.AgglomerativeClustering(n_clusters=8, linkage='ward')
complete = cluster.AgglomerativeClustering(n_clusters=8, linkage='complete')
average = cluster.AgglomerativeClustering(n_clusters=8, linkage='average')
single = cluster.AgglomerativeClustering(n_clusters=8, linkage='single')

clustering_algorithms = (
        ('Single Linkage', single),
        ('Average Linkage', average),
        ('Complete Linkage', complete),
        ('Ward Linkage', ward),
    )
clusterings = [single,average, complete,ward]

fig, axes = plt.subplots(2,2, figsize=(16,16))
for i, algorithm in enumerate(clusterings):
        y_pred = algorithm.fit_predict(data_pca)
        #print(i)
        ax = axes.flatten()[i]
        scprep.plot.scatter2d(data_phate, c=y_pred, cmap=cluster_cmap,
                          title='{} - ({})'.format(clustering_algorithms[i][0], len(np.unique(y_pred))), 
                          ticks=False, label_prefix="PHATE", legend=False, discrete=True,
                           ax=ax)

plt.savefig("hiecey.png")
plt.savefig("hiecey.pdf")

plt.figure(figsize=(10,10))
scprep.plot.scatter2d(data_phate, c=kmeans_clusters, cmap=cluster_cmap,
            title='{} - ({})'.format("Kmeans", len(np.unique(kmeans_clusters))), 
            ticks=False, label_prefix="PHATE", legend=False, discrete=True)
plt.savefig("kmeans.pdf")

for i, algorithm in enumerate(clustering_algorithms):
    ax = axes.flatten()[i]
    scprep.plot.scatter2d(data_phate, c=clusterings[algorithm], cmap=cluster_cmap,
                          title='{} - ({})'.format(algorithm, len(np.unique(clusterings[alg]))), 
                      ticks=False, label_prefix="PHATE", legend=False, discrete=True,
                         ax=ax)

os.chdir("/mnt/d/SingleCell/Workshop/workshop/02.Clustering/data/")
data_pca = scprep.reduce.pca(adata.X, n_components=40)
cluster_cmap = plt.cm.tab20(np.linspace(0, 1, 20))
phate_op = phate.PHATE(knn=5, n_jobs=-2)
data_phate = phate_op.fit_transform(data_pca)

fig, axes = plt.subplots(2,3, figsize=(16,16))
clu=[5,6,7,8,9,10]
for i,k in enumerate(clu):
    ax = axes.flatten()[i]
    phenograph_clusters, _, _ = phenograph.cluster(data_pca,k=k)#default
    scprep.plot.scatter2d(data_phate, c=phenograph_clusters, cmap=cluster_cmap,
            title='{} - ({})-(k={})'.format("phenograph", len(np.unique(phenograph_clusters)),k), 
            ticks=False, label_prefix="PHATE", legend=False, discrete=True,ax=ax)


plt.savefig("phenograpall.png")

G = gt.Graph(scprep.reduce.pca(data, n_components=50))
G_igraph = G.to_igraph()

tic = time.time()
part = louvain.find_partition(G_igraph, louvain.RBConfigurationVertexPartition, 
                              weights="weight", resolution_parameter=1)

louvain_clusters = np.array(part.membership)
print('Finished Louvain in {:.2f} seconds'.format(time.time() - tic))

spec_op = sklearn.cluster.SpectralClustering(n_clusters=20, affinity='precomputed')
spectral_clusters = spec_op.fit_predict(G.K.toarray())
print('Finished Spectral clustering in {:.2f} seconds'.format(time.time() - tic))

clusterings = {'Phenograph':phenograph_clusters,
               'Louvain':louvain_clusters, 
               'KMeans':kmeans_clusters, 
               'Spectral':spectral_clusters}

for alg in clusterings:
    cl_nu = scprep.utils.sort_clusters_by_values(clusterings[alg], -data_phate.iloc[:,0])
    clusterings[alg] = cl_nu

cluster_cmap = plt.cm.tab20(np.linspace(0, 1, 20))


fig, axes = plt.subplots(2,2, figsize=(16,16))

for i, algorithm in enumerate(clusterings):
    ax = axes.flatten()[i]
    scprep.plot.scatter2d(data_phate, c=clusterings[algorithm], cmap=cluster_cmap,
                          title='{} - ({})'.format(algorithm, len(np.unique(clusterings[alg]))), 
                      ticks=False, label_prefix="PHATE", legend=False, discrete=True,
                         ax=ax)

plt.savefig("clustering.pdf")

n_rows = 10
n_cols = 3

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*4,n_rows*4))
axes = axes.flatten()
clusters = clusterings['Spectral']
for i, ax in enumerate(axes):
    try:
        curr_cluster = np.unique(clusters)[i]
    except IndexError:
        ax.axis('off')
        continue
    # Returns array([False, True,...,False]) indicating if each cell is in the
    # current cluster
    curr_mask = clusters == curr_cluster  
    scprep.plot.scatter2d(data_phate.loc[~curr_mask] , color='grey', zorder=0, s=1, ax=ax)
    scprep.plot.scatter2d(data_phate.loc[curr_mask], color=cluster_cmap[curr_cluster], title='Cluster {}'.format(curr_cluster),
                          ticks=False, label_prefix='PHATE ', ax=ax)


fig.tight_layout()

plt.savefig("all.clustering.pdf")