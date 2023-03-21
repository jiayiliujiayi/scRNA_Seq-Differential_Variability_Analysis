setwd("~/denSNE/10.0.illustration/realCOVID/01.preprocessing")
library(Seurat)
source("~/pkg/densne/densne.R")

library(Seurat)
library(magrittr)
library(dplyr)
library(data.table)
library(edgeR)
library(RColorBrewer)
library(stringr)
library(cowplot)
library(ggplot2)
library(future)
library(parallel)
library(ggpubr)

datetag = gsub("-", "", Sys.Date())

source("~/denSNE/8.0calculate-affinity_batches/01normalization.R")
source("~/denSNE/8.0calculate-affinity_batches/02centering_scaling.R")
source("~/denSNE/8.0calculate-affinity_batches/03PCA.R")

# import raw
counts = Matrix::readMM('../99.processed/counts_MildSevere.mtx')
meta = fread("../99.processed/metadata_MildSevere.csv")
genes = read.delim("../99.processed/features.csv", header = 1, sep = ",")

# data cleaning
## init Seurat
### counts and metadata
counts.mat = as.matrix(counts)
rownames(counts.mat) = meta$covid_index
colnames(counts.mat) = genes$X
subgenes = genes %>% filter(feature_types == "Gene Expression")
counts.mat = counts.mat[, subgenes$X]
counts.mat = t(counts.mat)

meta$cell = meta$covid_index

all.equal(colnames(counts.mat), meta$cell)
rownames(meta) = meta$cell

### obj
obj = CreateSeuratObject(counts.mat, meta.data = meta)
###QC
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
obj$cell_type = obj$initial_clustering

obj@meta.data %>% select(cell_type, Status_on_day_collection_summary) %>% table

# write out
saveRDS(obj, "../99.processed/all_MildSevere.obj")
