.libPaths("/home/XXX/R/rocker-rstudio/4.0/")

source("~/pkg/densne/densne.R")
library(future)
library(parallel)
library(GDAtools)
library(rdist)
library(dplyr)     
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggplot2)
library(Matrix, lib.loc = "/home/XXX/R/rocker-rstudio/4.0/")
dyn.load('/home/XXX/miniconda3/lib/libxml2.so.2'); 
dyn.load('/home/XXX/R/rocker-rstudio/4.0/igraph/libs/igraph.so');
dyn.load('/home/XXX/miniconda3/lib/libglpk.so.40');library(Seurat)

library(magrittr)
library(igraph)

datetag = gsub("\\-", "", Sys.Date())

# extract args
args = commandArgs(trailingOnly=T)
minfolder = input_datetag = args[grep('--minfolder=', args)] %>% gsub("--minfolder=", "", .)
maxfolder = input_datetag = args[grep('--maxfolder=', args)] %>% gsub("--maxfolder=", "", .)
input_datetag = args[grep('--input_datetag=', args)] %>% gsub("--input_datetag=", "", .) %>% as.character
gene_perc = args[grep('--gene_perc=', args)] %>% gsub("--gene_perc=", "", .) %>% as.numeric
factor_upper = args[grep('--factor_upper=', args)] %>% gsub("--factor_upper=", "", .) %>% as.numeric
factor_lower = args[grep('--factor_lower=', args)] %>% gsub("--factor_lower=", "", .) %>% as.numeric
ncells_per_type = args[grep('--ncells_per_type=', args)] %>% gsub("--ncells_per_type=", "", .) %>% as.numeric
# HERE Change n Types
ncells = 6 * ncells_per_type

path.in = sprintf("/home/XXX/denSNE/8.0calculate-affinity_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell", 
                  input_datetag, gene_perc, 
                  factor_lower,
                  factor_upper, 
                  ncells)
files = list.files(path.in, pattern = "method15.txt", recursive = T)
filesnumber = files %>% gsub("/method15.txt",  "", .) %>% as.numeric()
not.avail = c(minfolder:maxfolder) %>% .[!. %in% filesnumber]
print(sprintf("not avail in %s", path.in))
print(not.avail)


mclapply(not.avail, mc.cores = 5,
         function(X){
           # source functions -------------------------------------------------
           setwd("~/denSNE/8.0calculate-affinity_batches/")
           source("../scripts/00importing.R")
           source("../scripts/01normalization.R")
           source("../scripts/02centering_scaling.R")
           source("../scripts/03PCA.R")
           source("../scripts/04calculating.affinities.R")
           # importing--------------------------------------------------
           path = sprintf("/home/XXX/denSNE/6.1simulation_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                          input_datetag, gene_perc, 
                          factor_lower,
                          factor_upper, 
                          ncells, X)
           print(sprintf("Start calculating for dataset %1.0f of %scells/type at %s", 
                         X, 
                         ncells, 
                         as.character(Sys.time())))
           # print(path)
           # print("importing datasets")
           imported = importData(path)
           obj = imported$obj; metadata = imported$metadata; counts.raw = imported$counts.raw
           
           # normalizing. centering and scaling -------------------------------------------
           # print("Normalizing")
           obj = norm_TP10K(obj)
           obj = scaleData(obj, cell.type = NULL)
           
           obj = FindVariableFeatures(obj, nfeatures = 2000)
           
           # PCA ---------------------------------------------------------------------
           # print("Performing PCA")
           obj = RunPCA(obj, assay = "RNA", features = VariableFeatures(obj), 
                        npcs = min(200, ncol(obj) - 1), 
                        verbose = F)
         
           # calculating affinities --------------------------------------------------
           # print("Calculating distances to medoids")
           affinity = calc_affinity(obj, method = "dist.medoid", cell.type = "Major.cell.type")
           colnames(affinity)[1] = "method15"
           
           dir.create(
             sprintf("/scratch/XXX/denSNE/8.0calculate-affinity_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell", 
                     input_datetag, gene_perc, factor_lower, factor_upper, 
                     ncells)
           )
           
           dir.create(
             sprintf("/scratch/XXX/denSNE/8.0calculate-affinity_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                     input_datetag, gene_perc, factor_lower, factor_upper, 
                     ncells, X)
           )
           
           outpath = sprintf("/scratch/XXX/denSNE/8.0calculate-affinity_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                             input_datetag, gene_perc, factor_lower, factor_upper, 
                             ncells, X)
           
           # print("Saving files")
           setwd(outpath)
           write.table(affinity, 
                       "method15.txt", 
                       quote = F, col.names = T, row.names = F, sep = "\t")
           
           print(sprintf("saved distances for dataset %1.0f of %scells/type, at %s", 
                         X, 
                         ncells, 
                         Sys.time() %>% as.character))
         })

# import ------------------------------------------------------------------
# #obj = readRDS("../6.1simulation_batches/out/20220125-var50_1.5_3.0-6type-600cell/1/set1.rds")
# #obj = CreateSeuratObject(obj)
# print("importing datasets")
# source("../scripts/00importing.R")
# 
# # normalizing -------------------------------------------------------------
# print("Normalizing")
# obj = norm_TP10K(obj)
# 
# # centering and scaling ---------------------------------------------------
# obj = scaleData(obj, cell.type = NULL)
# 
# # PCA ---------------------------------------------------------------------
# ## find variable genes: TP10K
# var.genes_tp = FindVariableFeatures(obj, nfeatures = 2000) %>% VariableFeatures()
# ## find variable genes: SCT
# obj.sct = Seurat::SCTransform(obj, variable.features.n = 2000, verbose = T)
# var.genes_sct = VariableFeatures(obj.sct)
# ## get variable genes 
# var.genes = intersect(var.genes_tp, var.genes_sct)
# 
# print("Performing PCA")
# obj = RunPCA(obj, assay = "RNA", features = var.genes, 
#              npcs = min(200, ncol(obj) - 1), 
#              verbose = F)
# 
# # calculating affinities --------------------------------------------------
# print("Calculating distances to medoids")
# affinity = calc_affinity(obj, method = "local.radius")
# colnames(affinity)[1] = "method15"
# 
# dir.create(
#   sprintf("./out/%s-var%.0f_%.1f_%.1f-6type-%scell", 
#           input_datetag, gene_perc, factor_lower, factor_upper, 
#           ncells)
# )
# 
# dir.create(
#   sprintf("./out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
#           input_datetag, gene_perc, factor_lower, factor_upper, 
#           ncells, folder)
# )
# 
# outpath = sprintf("./out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
#                   input_datetag, gene_perc, factor_lower, factor_upper, 
#                   ncells, folder)
# 
# print("Saving files")
# write.table(affinity, 
#             sprintf("%s/method15.txt", outpath), 
#             quote = F, col.names = T, row.names = F, sep = "\t")
# 
# print("finished")