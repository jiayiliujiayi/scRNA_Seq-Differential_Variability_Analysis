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
dyn.load('~/.conda/envs/jlenv/lib/libxml2.so.2'); 
dyn.load('/home/XXX/R/rocker-rstudio/4.0/igraph/libs/igraph.so');
dyn.load('~/.conda/envs/jlenv/lib/libglpk.so.40');library(Seurat)

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
files = list.files(path.in, pattern = "method41.txt", recursive = T)
filesnumber = files %>% gsub("/method41.txt",  "", .) %>% as.numeric()
not.avail = c(minfolder:maxfolder) %>% .[!. %in% filesnumber]
print(sprintf("not avail in %s", path.in))
print(not.avail)

lapply(not.avail, #mc.cores = 5,
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
           
           
           # normalizing -------------------------------------------------------------
           # print("Normalizing")
           obj = NormalizeData(obj, verbose = F)
           
           # centering and scaling ---------------------------------------------------
           obj = ScaleData(obj, verbose = F)
           
           # PCA ---------------------------------------------------------------------
           # print("Performing PCA")
           obj = do_PCA(obj, 
                        PCA.input = "all.genes")
           
           # denSNE ---------------------------------------------------------------------
           coord_pca = Embeddings(obj, reduction = "pca")
           coord_pca = coord_pca[, 1:min(200, ncol(coord_pca))]
           out_densne = run_densne(coord_pca, no_dims=2, perplexity=50,
                                   theta=0.5,
                                   randseed=X,
                                   verbose=FALSE, use_pca=FALSE,
                                   max_iter=1000, dens_frac=0.3, # as in the paper
                                   final_dens=TRUE)
           densne_embeddings = out_densne[[1]]
           colnames(densne_embeddings) = paste0("denSNE_", 1:2)
           rownames(densne_embeddings) = colnames(obj)
           obj[["denSNE"]] <- CreateDimReducObject(embeddings = densne_embeddings,
                                                   key = "denSNE_", 
                                                   assay = DefaultAssay(obj))
           
           # calculating affinities --------------------------------------------------
           # print("Calculating distances to medoids")
           affinity = calc_affinity(obj, method = "2d.denSNE_to_DM", cell.type = "Major.cell.type")
           colnames(affinity)[1] = "method41"
           
           
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
                       "method41.txt", 
                       quote = F, col.names = T, row.names = F, sep = "\t")
           
           # print("finished")           
           print(sprintf("saved distances for dataset %1.0f of %scells/type, at %s", 
                         X, 
                         ncells, 
                         Sys.time() %>% as.character))
           
           # remove previous session outputs
           folder <- tempdir()
           # get all files in the directories, recursively
           f <- list.files(folder, include.dirs = F, full.names = T, recursive = T)
           # remove the files
           file.remove(f)
         })

# # normalizing. centering and scaling --------------------------------------
# print("Normalizing")
# obj = NormalizeData(obj)
# 
# # centering and scaling ---------------------------------------------------
# obj = ScaleData(obj)
# 
# # find variable features --------------------------------------------------
# obj = FindVariableFeatures(obj)
# 
# # PCA ---------------------------------------------------------------------
# print("Performing PCA")
# obj = RunPCA(obj, npcs = min(200, ncol(obj) - 1))
# 
# # calculating affinities --------------------------------------------------
# print("Calculating distances to medoids")
# affinity = calc_affinity(obj, method = "local.radius")
# colnames(affinity)[1] = "method41"
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
#             sprintf("%s/method41.txt", outpath), 
#             quote = F, col.names = T, row.names = F, sep = "\t")
# 
# print("finished")
