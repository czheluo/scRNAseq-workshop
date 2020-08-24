#install.packages("BiocManager")
#BiocManager::install("slingshot")

library(slingshot)
library(SingleCellExperiment)
#install.packages("devtools")
#devtools::install_github("dynverse/dyno")
#devtools::install_github("dynverse/dyno", host = "http://api.github.com")
library(babelwhale)
config <- create_docker_config()
set_default_config(config)
library(dyno)
library(tidyverse)
#setwd("I:\\MAJORBIO\\MAJORBIO\\SingleCell\\SingleCellWorkshop")
data("fibroblast_reprogramming_treutlein")
#load("C:\\Users\\Administrator\\Desktop\\dyno-master\\data/fibroblast_reprogramming_treutlein.rda")
dataset <- wrap_expression(
  counts = fibroblast_reprogramming_treutlein$counts,
  expression = fibroblast_reprogramming_treutlein$expression
)

guidelines <- guidelines(
  dataset,
  answers = answer_questions(
    dataset,
    multiple_disconnected = FALSE,
    expect_topology = TRUE,
    expected_topology = "bifurcation"
  )
)
answers <- dynguidelines::answer_questions(
  multiple_disconnected = NULL, 
  expect_topology = NULL, 
  expected_topology = NULL, 
  n_cells = 392, 
  n_features = 2000, 
  memory = "2GB", 
  prior_information = c("start_id", "end_id", "end_n", "start_n", "leaves_n", "groups_n", "features_id", "dimred"), 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)
guideliness <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected
model_slingshot <- infer_trajectory(dataset, methods_selected[1])
model_paga_tree <- infer_trajectory(dataset, methods_selected[2])
model_paga <- infer_trajectory(dataset, methods_selected[4])
model_grandprix <- infer_trajectory(dataset, methods_selected[3])
model <- model_slingshot

plot_dimred(
  model_slingshot,
  expression_source = dataset$expression,
  grouping = dataset$grouping
)

add_prior_information(dataset, end_n = 1)

model <- infer_trajectory(dataset, "slingshot")
model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
write.csv(fibroblast_reprogramming_treutlein$grouping, file='Treutlein.groups.csv')
write.csv(as.matrix(dataset$counts), file='Treutlein.counts.csv')

write.csv(as.matrix(dataset$expression), file='Treutlein.expression.csv')
write.csv(as.matrix(dataset$feature_info), file='Treutlein.gene_names.csv')

# Load counts, PHATE, and cluster labels
counts <- as.data.frame(read.csv('Treutlein.counts.csv'))
phate = as.matrix(read.csv('Treutlein.expression.csv'))
pclusters = read.csv('Treutlein.groups.csv')

# Create SingleCellExperiment
# How am I actually supposed to have the column names not be passed as genes??
sim <- SingleCellExperiment(assays = List(counts = t(as.matrix(counts[,2:2001]))))

# Add dim red data and clusters to SCE
reducedDims(sim) <- SimpleList(PHATE=phate)
colData(sim)$pclusters <- pclusters[,2]

# Optional, plot data and clusters
library(RColorBrewer)
plot(phate, col = brewer.pal(5,"Set1")[sim$pclusters], pch=16, asp = 1)

# Do Slingshot
sce <- slingshot(sim, clusterLabels = 'pclusters', reducedDim = 'PHATE')

summary(sce$slingPseudotime_3)

# Plot Slingshot
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

# For some reason not all points are plotted here
plot(reducedDims(sce)$PHATE, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

# This plots the 'scaffold'
plot(reducedDims(sce)$PHATE, col = brewer.pal(5,"Set1")[sim$pclusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

#You can get the orderings of the points from these variables
sce$slingPseudotime_1
sce$slingPseudotime_2

model <- model_paga
plot_dimred(
  model,
  expression_source = dataset$expression,
  grouping = dataset$grouping
)
