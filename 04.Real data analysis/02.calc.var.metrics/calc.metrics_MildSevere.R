setwd("xxxx")
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

source("../02.Preprocessing & Calculating variability metrics/scripts/01normalization.R")
source("../02.Preprocessing & Calculating variability metrics/scripts/02centering_scaling.R")
source("../02.Preprocessing & Calculating variability metrics/scripts/03PCA.R")
source("../02.Preprocessing & Calculating variability metrics/scripts/04calculating.affinities.R")

# import obj.immune list -----
x <- readRDS("../99.processed/all_MildSevere.obj")
x$Major.cell.type = x$initial_clustering
x$sim.method = x$Status_on_day_collection_summary
x$celltype_status = paste0(x$Major.cell.type, "_", x$Status_on_day_collection_summary)

# subset
nCell_per.type = table(x$celltype_status)
index = vector(mode = "list", length = length(nCell_per.type)) %>% `names<-`(names(nCell_per.type))
for (i in 1:length(index)) {
  celltype = names(index)[i]
  if (nCell_per.type[celltype] <= 5000) {
    index[[celltype]] = which(x$celltype_status == celltype)
  } else{
    set.seed(41); index[[celltype]] = which(x$celltype_status == celltype) %>% .[sample(1:length(.), 000)]
  } #else {}
}
index = unlist(index) %>% unname
obj = x[, index]
#cellid = colnames(obj)[index]; cell_type = obj$cell.type[index] #cell_type %>% table

# calc metrics -----
## TP_A_LR/DM -----
obj = norm_TP10K(obj)#NormalizeData(obj, normalization.method = "RC", scale.factor = 1000)
obj = scaleData(obj)#ScaleData(x, cell.type = NULL)
obj = do_PCA(obj, PCA.input = "all.genes")

saveRDS(obj, "../99.processed/MildSevere_withPCA.RDS")
obj = readRDS("../99.processed/MildSevere_withPCA.RDS")
### denSNE ----
coord_pca = Embeddings(obj, reduction = "pca")
coord_pca = coord_pca[, 1:min(200, ncol(coord_pca))]
out_densne = run_densne(coord_pca, no_dims=2, perplexity=50,
                        theta=0.5,
                        randseed=-1,
                        verbose=FALSE, use_pca=FALSE,
                        max_iter=1000, dens_frac=0.3, # as in the paper
                        final_dens=TRUE)
densne_embeddings = out_densne[[1]]
colnames(densne_embeddings) = paste0("denSNE_", 1:2)
rownames(densne_embeddings) = colnames(obj)
obj[["denSNE"]] <- CreateDimReducObject(embeddings = densne_embeddings,
                                        key = "denSNE_", 
                                        assay = DefaultAssay(obj))

TP_A_denSNE = calc_affinity(obj, method = "2d.denSNE_to_DM", cell.type = "Major.cell.type")

saveRDS(TP_A_denSNE, "./MildSevere_TP_A_denSNE.RDS")

