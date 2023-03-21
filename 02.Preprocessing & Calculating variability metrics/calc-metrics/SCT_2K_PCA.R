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
files = list.files(path.in, pattern = "method14.txt", recursive = T)
filesnumber = files %>% gsub("/method14.txt",  "", .) %>% as.numeric()
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
  obj = SCTransform(obj, variable.features.n = 2000, verbose = F)
  
  # PCA ---------------------------------------------------------------------
  # print("Performing PCA")
  obj = RunPCA(obj, assay = "SCT", features = VariableFeatures(obj), 
               npcs = min(200, ncol(obj) - 1), 
               verbose = F)
  
  # calculating affinities --------------------------------------------------
  # print("Calculating distances to medoids")
  affinity = calc_affinity(obj, method = "dist.medoid", cell.type = "Major.cell.type")
  colnames(affinity)[1] = "method14"
  
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
  write.table(affinity, 
              sprintf("%s/method14.txt", outpath), 
              quote = F, col.names = T, row.names = F, sep = "\t")
  
  print(sprintf("saved distances for dataset %1.0f of %scells/type, at %s", 
                X, 
                ncells, 
                 
                Sys.time() %>% as.character))
         })
