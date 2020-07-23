import scanpy.api as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from maren_codes import combat as c
from gprofiler import gprofiler

import warnings
from rpy2.rinterface import RRuntimeWarning
from rpy2.robjects import pandas2ri

# Perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain_r1')
sc.tl.louvain(adata, resolution=0.5, key_added='louvain_r0.5')



### scripe
import sys
import os
import time
#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
#font = {'size'   : 14}
#mpl.rc('font', **font)
import numpy as np
import pandas as pd
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

data = pd.read_pickle('embryoid_body_data.pickle.gz')
metadata = pd.DataFrame([ix.split('_')[1] for ix in data.index], columns=['sample'], index=data.index)
full_marker_genes = ['ARID3A (ENSG00000116017)', 'ASCL2 (ENSG00000183734)',  'CD34 (ENSG00000174059)',
 'CDX2 (ENSG00000165556)', 'CER1 (ENSG00000147869)', 'DLX1 (ENSG00000144355)',
 'DMRT3 (ENSG00000064218)', 'EN2 (ENSG00000164778)', 'EOMES (ENSG00000163508)',
 'FOXA2 (ENSG00000125798)', 'FOXD3-AS1 (ENSG00000230798)', 'GATA3-AS1 (ENSG00000197308)',
 'GATA4 (ENSG00000136574)', 'GATA5 (ENSG00000130700)', 'GATA6-AS1 (ENSG00000266010)',
 'GBX2 (ENSG00000168505)', 'GLI3 (ENSG00000106571)', 'HOXA2 (ENSG00000105996)',
 'HOXB1 (ENSG00000120094)', 'HOXB4 (ENSG00000182742)', 'HOXD13 (ENSG00000128714)',
 'HOXD9 (ENSG00000128709)', 'ISL1 (ENSG00000016082)', 'KLF5 (ENSG00000102554)',
 'KLF7 (ENSG00000118263)', 'LEF1 (ENSG00000138795)', 'LHX2 (ENSG00000106689)',
 'LHX5 (ENSG00000089116)', 'LMX1A (ENSG00000162761)', 'MAP2 (ENSG00000078018)',
 'MIXL1 (ENSG00000185155)', 'MYCBP (ENSG00000214114)', 'NANOG (ENSG00000111704)',
 'NES (ENSG00000132688)', 'NKX2-1 (ENSG00000136352)', 'NKX2-5 (ENSG00000183072)',
 'NKX2-8 (ENSG00000136327)', 'NPAS1 (ENSG00000130751)', 'NR2F1-AS1 (ENSG00000237187)',
 'OLIG1 (ENSG00000184221)', 'OLIG3 (ENSG00000177468)', 'ONECUT1 (ENSG00000169856)',
 'ONECUT2 (ENSG00000119547)', 'OTX2 (ENSG00000165588)', 'PAX3 (ENSG00000135903)',
 'PAX6 (ENSG00000007372)', 'PDGFRA (ENSG00000134853)', 'PECAM1 (ENSG00000261371)',
 'POU5F1 (ENSG00000204531)', 'SATB1 (ENSG00000182568)', 'SIX2 (ENSG00000170577)',
 'SIX3-AS1 (ENSG00000236502)', 'SIX6 (ENSG00000184302)', 'SOX13 (ENSG00000143842)',
 'SOX10 (ENSG00000100146)', 'SOX15 (ENSG00000129194)', 'SOX17 (ENSG00000164736)',
 'SOX9 (ENSG00000125398)', 'TTLL10 (ENSG00000162571)', 'TAL1 (ENSG00000162367)',
 'TBX15 (ENSG00000092607)', 'TBX18 (ENSG00000112837)', 'TBX5 (ENSG00000089225)',
 'TNNT2 (ENSG00000118194)', 'WT1 (ENSG00000184937)', 'ZBTB16 (ENSG00000109906)',
 'ZIC2 (ENSG00000043355)', 'ZIC5 (ENSG00000139800)', 'ACTB (ENSG00000075624)',
 'HAND1 (ENSG00000113196)']
import magic
data_magic = magic.MAGIC().fit_transform(data, genes=full_marker_genes)

data_phate = phate.PHATE().fit_transform(data)
# alternative: umap.UMAP(), sklearn.manifold.TSNE()
data_phate = pd.DataFrame(data_phate, index=data.index)
plt.figure(figsize=(10,10))
scprep.plot.scatter2d(data_phate, c=metadata['sample'], figsize=(12,8), cmap="Spectral",
                      ticks=False, label_prefix="PHATE")
plt.savefig("phatedata.pdf")
home = os.path.expanduser('./')
file_path = os.path.join(home, 'EBT_counts.pkl.gz')
if not os.path.exists(file_path):
    scprep.io.download.download_google_drive(id='1Xz0ONnRWp2MLC_R6r74MzNwaZ4DkQPcM',
                        destination=os.path.dirname(file_path))

data = pd.read_pickle(file_path)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline

import sklearn
import sklearn.cluster
import sklearn.manifold

import os
import tasklogger
import phate
import umap

import graphtools as gt
import magic
import phenograph
import louvain
data_pca = scprep.reduce.pca(data, n_components=50, method='dense')

phenograph_clusters, _, _ = phenograph.cluster(data_pca)

phenograph_clusters, _, _ = phenograph.cluster(data_pca1)

with tasklogger.log_task("KMeans"):
    kmeans_clusters = sklearn.cluster.KMeans(n_clusters=20).fit_predict(data_pca)

G = gt.Graph(data_pca)
G_igraph = G.to_igraph()

with tasklogger.log_task("Louvain"):
    partition = louvain.find_partition(G_igraph, louvain.RBConfigurationVertexPartition, 
                                       weights="weight", resolution_parameter=1)
    louvain_clusters = np.array(partition.membership)

with tasklogger.log_task("Spectral clustering"):
    spec_op = sklearn.cluster.SpectralClustering(n_clusters=20, affinity='precomputed')
    spectral_clusters = spec_op.fit_predict(G.K)

clusterings = {'Phenograph':phenograph_clusters,
               'Louvain':louvain_clusters, 
               'KMeans':kmeans_clusters, 
               'Spectral':spectral_clusters}

for alg in clusterings:
    cl_nu = scprep.utils.sort_clusters_by_values(clusterings[alg], -data_phate.iloc[:,0])
    clusterings[alg] = cl_nu

cluster_cmap = plt.cm.tab20(np.linspace(0, 1, 21))
fig, axes = plt.subplots(2,2, figsize=(16,16))

for i, algorithm in enumerate(clusterings):
    ax = axes.flatten()[i]
    scprep.plot.scatter2d(data_phate, c=clusterings[algorithm], cmap=cluster_cmap,
                          title='{} - ({})'.format(algorithm, len(np.unique(clusterings[alg]))), 
                      ticks=False, label_prefix="PHATE", legend=False, discrete=True,
                         ax=ax)

plt.savefig("fourallcluster.pdf")
plt.savefig("fourallcluster.png")

cluster_cmap = plt.cm.tab20(np.linspace(0, 1, 22))
n_rows = 8
n_cols = 3

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*4,n_rows*4))
axes = axes.flatten()
clusters = clusterings['Louvain']#Spectral

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
plt.savefig("all.pdf")
plt.savefig("all.png")

all_clusterings = []
all_algorithms = list(clusterings.keys())
for algo in all_algorithms:
    all_clusterings.append(clusterings[algo])    
all_clusterings = np.vstack(all_clusterings)
dists = scipy.spatial.distance.pdist(all_clusterings, metric=sklearn.metrics.adjusted_rand_score)
# squareform assumes diagonals will be 0, but they're actually 1 because this is a similarity metric
# so we need to add 1's on the diagonal with np.eye()
dists = scipy.spatial.distance.squareform(dists) + np.eye(4)
sns.clustermap(dists, xticklabels=all_algorithms, yticklabels=all_algorithms)

fig, ax = plt.subplots(1, figsize=(12, 5))

curr_gene = 'POU5F1'
curr_expression = scprep.select.select_cols(EBT_counts, starts_with=curr_gene)
curr_expression = curr_expression.iloc[:,0]
scprep.plot.jitter(clusters, curr_expression, c=clusters, cmap=cluster_cmap, ax=ax,
                  legend_anchor=(1,1))
ax.set_title(curr_gene)


fig, ax = plt.subplots(1, figsize=(12, 5))

curr_gene = 'NANOG'
curr_expression = scprep.select.select_cols(EBT_counts, starts_with=curr_gene)
curr_expression = curr_expression.iloc[:,0]
scprep.plot.jitter(clusters, curr_expression, c=clusters, cmap=cluster_cmap, ax=ax,
                  legend_anchor=(1,1))
ax.set_title(curr_gene)

##geneplot
fig, axes = plt.subplots(1,3, figsize=(14,4))
axes = axes.flatten()


genes_for_plotting = ['NANOG', 'POU5F1', 'HAND1']

for i, ax in enumerate(axes):
    curr_gene = genes_for_plotting[i]
    

    expression = scprep.select.select_cols(data, starts_with=curr_gene).to_dense()
    if expression.shape[1] > 1:
        expression = expression[expression.columns[0]]
        expression = pd.DataFrame(expression)
    
    sort_index = expression.sort_values(by=expression.columns[0]).index
    
    scprep.plot.scatter2d(data_phate.loc[sort_index], c=expression.loc[sort_index], shuffle=False,
                         title=curr_gene, ticks=None, label_prefix='PHATE ', ax=ax)
    
fig.tight_layout()


