
setwd("D:\\SingleCell\\Workshop\\workshop\\01.QC\\session-qc\\data")

suppressMessages(require(Seurat))
suppressMessages(require(scater))
suppressMessages(require(Matrix))
v3.1k <- Read10X_h5("pbmc_1k_v3_filtered_feature_bc_matrix.h5", use.names = T)
v2.1k <- Read10X_h5("pbmc_1k_v2_filtered_feature_bc_matrix.h5", use.names = T)
p3.1k <- Read10X_h5("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5", use.names = T)

## Genome matrix has multiple modalities, returning a list of matrices for this genome
# select only gene expression data from the CITE-seq data.
p3.1k <- p3.1k$`Gene Expression`


#First, create Seurat objects for each of the datasets, and then merge into one large seurat object.
sdata.v2.1k <- CreateSeuratObject(v2.1k, project = "v2.1k")
sdata.v3.1k <- CreateSeuratObject(v3.1k, project = "v3.1k")
sdata.p3.1k <- CreateSeuratObject(p3.1k, project = "p3.1k")

# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the datasets.
alldata <- merge(sdata.v2.1k, c(sdata.v3.1k,sdata.p3.1k), add.cell.ids=c("v2.1k","v3.1k","p3.1k"))

# also add in a metadata column that indicates v2 vs v3 chemistry
chemistry <- rep("v3",ncol(alldata))
chemistry[Idents(alldata) == "v2.1k"] <- "v2"
alldata <- AddMetaData(alldata, chemistry, col.name = "Chemistry")
alldata

## An object of class Seurat
## 33538 features across 2931 samples within 1 assay
## Active assay: RNA (33538 features)
# check number of cells from each sample, is stored in the orig.ident slot of metadata and is autmatically set as active ident.
table(Idents(alldata))

##
## p3.1k v2.1k v3.1k
##   713   996  1222


##Calculate mitochondrial proportion

#Seurat automatically calculates some QC-stats,
#like number of UMIs and features per cell. Stored in columns nCount_RNA & nFeature_RNA of the metadata.

head(alldata@meta.data)

#We will manually calculate the proportion of mitochondrial reads and add to the metadata table.

mt.genes <- rownames(alldata)[grep("^MT-",rownames(alldata))]
C<-GetAssayData(object = alldata, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
alldata <- AddMetaData(alldata, percent.mito, col.name = "percent.mito")

###Calculate ribosomal proportion

#In the same manner we will calculate the proportion gene expression that comes from ribosomal proteins. NOTE - add text on why!

rb.genes <- rownames(alldata)[grep("^RP[SL]",rownames(alldata))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
alldata <- AddMetaData(alldata, percent.ribo, col.name = "percent.ribo")



#plotQC
VlnPlot(alldata, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(alldata, features = "nCount_RNA", pt.size = 0.1) + NoLegend()

VlnPlot(alldata, features = "percent.mito", pt.size = 0.1) + NoLegend()

VlnPlot(alldata, features = "percent.ribo", pt.size = 0.1) + NoLegend()


#As you can see, the v2 chemistry gives lower gene detection, but higher detection of ribosomal proteins. As the ribosomal proteins are highly expressed they will make up a larger proportion of the transcriptional landscape when fewer of the lowly expressed genes are detected.

#And we can plot the different QC-measures as scatter plots


FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(alldata, feature1 = "nFeature_RNA", feature2 = "percent.mito")

FeatureScatter(alldata, feature1="percent.ribo", feature2="nFeature_RNA")


p<-FeatureScatter(alldata, feature1="percent.ribo", feature2="percent.mito",
                  cells = WhichCells(alldata, expression = orig.ident == "v3.1k"))
FeatureScatter(alldata, feature1="percent.ribo", feature2="percent.mito",
               cells = WhichCells(alldata, expression = orig.ident == "v2.1k"))
ggsave(p,filename ="ribo.mito.pdf" )
#We can also subset the data to only plot one sample.

FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               cells = WhichCells(alldata, expression = orig.ident == "v2.1k") )


#Filtering
#Mitochondrial filtering

#select cells with percent.mito < 25
selected <- WhichCells(alldata, expression = percent.mito < 25)
length(selected)
# and subset the object to only keep those cells
data.filt <- subset(alldata, cells = selected)
data.filt <- subset(alldata, subset = percent.mito < 25)
# plot violins for new data
VlnPlot(data.filt, features = "percent.mito")


#Gene detection filtering


#start with cells with many genes detected.
high.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA > 4100)
high.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA > 2000 & orig.ident == "v2.1k")

# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(high.det.v2,high.det.v3)))

# check number of cells
ncol(data.filt)

#Filter the cells with low gene detection (low quality libraries) with less than 1000 genes for v2 and < 500 for v2.

#start with cells with many genes detected.
low.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA < 1000 & orig.ident != "v2.1k")
low.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA < 500 & orig.ident == "v2.1k")

# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.det.v2,low.det.v3)))

# check number of cells
ncol(data.filt)


#Plot QC-stats again
VlnPlot(data.filt, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()

VlnPlot(data.filt, features = "nCount_RNA", pt.size = 0.1) + NoLegend()

VlnPlot(data.filt, features = "percent.mito", pt.size = 0.1) + NoLegend()

VlnPlot(data.filt, features = "percent.ribo", pt.size = 0.1) + NoLegend()

# and check the number of cells per sample before and after filtering
table(Idents(alldata))

table(Idents(data.filt))


#Calculate cell-cycle scores


#Seurat has a function for calculating cell cycle scores based on a list of know S-phase and G2/M-phase genes.

data.filt <- CellCycleScoring(
  object = data.filt,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

VlnPlot(data.filt, features = c("S.Score","G2M.Score"))


#Scater

sce <- as.SingleCellExperiment(data.filt)


#Calculate QC-metrics

# calculate all qc-metrics
sce <- calculateQCMetrics(sce, feature_controls = list(mito = mt.genes))

# check what all entries are -
colnames(colData(sce))

colnames(rowData(sce))

#Most expressed features
plotHighestExprs(sce, exprs_values = "counts")


#Cumulative expression

# plot each sample separately
plotScater(sce, block1 = "ident", nfeatures = 1000)

#Plot gene stats
plotRowData(sce, x = "n_cells_by_counts", y = "mean_counts")

#Plot cell stats
#In the same manner plotColData can plot any of the qc-measures for cells.

p1 <- plotColData(sce, x = "total_counts",
    y = "total_features_by_counts", colour_by = "ident")
p2 <- plotColData(sce, x = "pct_counts_feature_control",
    y = "total_features_by_counts", colour_by = "ident")
p3 <- plotColData(sce, x = "pct_counts_feature_control",
    y = "pct_counts_in_top_50_features", colour_by = "ident")
multiplot(p1, p2, p3, cols = 2)


#Identify outliers in QC-stats
#On method of identifying low quality cells is to run PCA on all the qc-stats and then identify outliers in PCA space.

sce <- runPCA(sce, use_coldata = TRUE,
    detect_outliers = TRUE)
## sROC 0.1-2 loaded

plotReducedDim(sce, use_dimred="PCA_coldata", colour_by = "ident")

# check if we have any outliers
table(colData(sce)$outlier)
plotPCA(sce,ncomponents=4,colour_by="ident")
plotPCA(sce,ncomponents=4,colour_by="percent.mito")

sce <- runUMAP(sce)
plotUMAP(object = sce, colour_by="ident")

p<-plotExplanatoryVariables(sce, variables =  c("ident","Chemistry","pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts", "S.Score","G2M.Score"))

data.filt <- CellCycleScoring(
  object = data.filt,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

VlnPlot(data.filt, features = c("S.Score","G2M.Score"))










