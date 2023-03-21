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
library(Matrix)
library(Seurat)
library(magrittr)
library(data.table)

datetag = gsub("\\-", "", Sys.Date())

# extract args
args = commandArgs(trailingOnly=T)
minfolder = args[grep('--minfolder=', args)] %>% gsub("--minfolder=", "", .)
maxfolder = args[grep('--maxfolder=', args)] %>% gsub("--maxfolder=", "", .)
input_datetag = args[grep('--input_datetag=', args)] %>% gsub("--input_datetag=", "", .) %>% as.character
gene_perc = args[grep('--gene_perc=', args)] %>% gsub("--gene_perc=", "", .) %>% as.numeric
factor_upper = args[grep('--factor_upper=', args)] %>% gsub("--factor_upper=", "", .) %>% as.numeric
factor_lower = args[grep('--factor_lower=', args)] %>% gsub("--factor_lower=", "", .) %>% as.numeric
ncells_per_type = args[grep('--ncells_per_type=', args)] %>% gsub("--ncells_per_type=", "", .) %>% as.numeric
# HERE Change n Types
ncells = 6 * ncells_per_type

path.in = sprintf("/home/XXXX/denSNE/8.2calculate-var_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell", 
                  input_datetag, gene_perc, 
                  factor_lower,
                  factor_upper, 
                  ncells)
#files1 = list.files(path.in, pattern = "MedianVar.txt", recursive = T)
files2 = list.files(path.in, pattern = "Median.Var.txt", recursive = T)
#files1number = files1 %>% gsub("/Median.Var.txt",  "", .) %>% as.numeric()
files2number = files2 %>% gsub("/Median.Var.txt",  "", .) %>% as.numeric()
filesnumber = c(files2number)
not.avail = c(minfolder:maxfolder) %>% .[!. %in% filesnumber]
print(sprintf("not avail in %s", path.in))
print(not.avail)


mclapply(not.avail, mc.cores = 7,
         function(X){
           # source functions -------------------------------------------------
           setwd("xxxx")
           source("xxxx/scripts/01normalization.R")
           # importing--------------------------------------------------
           path = sprintf("/home/XXXX/denSNE/6.1simulation_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                          input_datetag, gene_perc, 
                          factor_lower,
                          factor_upper, 
                          ncells, X)
           print(sprintf("Start calculating for dataset %1.0f of %scells/type at %s", 
                         X, 
                         ncells, 
                         as.character(Sys.time())))
           # importing ---------------------------------------------------------------
           files = list.files(path)
           copula1.path = grep("copula_result.set1", files)
           copula2.path = grep("copula_result.set2", files)
           copula.set1 = readRDS(sprintf("%s/%s", path, files[copula1.path[1]]))
           copula.set2 = readRDS(sprintf("%s/%s", path, files[copula2.path[1]]))
           counts = readRDS(sprintf("%s/combined.rds", path))
           
           # get altered gene id -----------------------------------------------------
           marg_param_into_df = function(X) {
             marg_param1 = data.frame(Gene = rownames(X$marginal_param1), 
                                      phi = X$marginal_param1[, 2], 
                                      mu = X$marginal_param1[, 3], 
                                      row.names = NULL)
             marg_param2 = data.frame(Gene = rownames(X$marginal_param2), 
                                      phi = X$marginal_param2[, 2], 
                                      mu = X$marginal_param2[, 3], 
                                      row.names = NULL)
             
             marg_param1$var = marg_param1$mu + marg_param1$mu^2 / marg_param1$phi
             marg_param2$var = marg_param2$mu + marg_param2$mu^2 / marg_param2$phi
             
             marg_param = rbind(marg_param1, marg_param2)
             marg_param = marg_param[!is.infinite(marg_param$phi), ]
             
             return(marg_param)
           }
           copula.set1_df.list = lapply(names(copula.set1), function(X) marg_param_into_df(copula.set1[[X]])) %>% `names<-`(names(copula.set1))
           copula.set1_df.list = lapply(names(copula.set1_df.list), function(X) mutate(copula.set1_df.list[[X]], CellType = gsub("\\..*", "", X))) %>% `names<-`(names(copula.set1))
           ## transform copula list into df for each cell type -------------------
           copula.set2_df.list = lapply(names(copula.set2), function(X) marg_param_into_df(copula.set2[[X]])) %>% `names<-`(names(copula.set2))
           copula.set2_df.list = lapply(names(copula.set2_df.list), function(X) mutate(copula.set2_df.list[[X]], CellType = gsub("\\..*", "", X))) %>% `names<-`(names(copula.set2))
           ## reduce copula list into one df -------------------
           copula.set1_df = Reduce(rbind, copula.set1_df.list)
           copula.set2_df = Reduce(rbind, copula.set2_df.list)
           
           var_copula.set1 = copula.set1_df %>% dplyr::select(., Gene, CellType, set1 = var)
           var_copula.set2 = copula.set2_df %>% dplyr::select(., Gene, CellType, set2 = var)
           
           var_copula.merged = Reduce(function(x, y) merge(x, y, by = c("Gene", "CellType")), 
                                      list(var_copula.set1, var_copula.set2))
           
           var_copula.merged$if.altered = ifelse(var_copula.merged$set1 == var_copula.merged$set2, "no", "yes")
           var_copula.merged$Gene.CellType = paste0(var_copula.merged$Gene, ".", var_copula.merged$CellType)
           ## get altered genes by cell type -------------------
           var_copula.merged_altered = var_copula.merged[var_copula.merged$if.altered == "yes", ]
           gene.celltype_altered = var_copula.merged_altered$Gene.CellType
           ## reform var_copula.merged_altered
           var_copula.merged_altered_reform = 
             var_copula.merged_altered %>% select(CellType, 
                                                  CD1 = set1, 
                                                  CD2 = set2, 
                                                  Gene.CellType)
           
           # calc val ------------------------------------------------------------
           ## normalizing -------------------------------------------------------------
           # print("Normalizing")
           obj = CreateSeuratObject(counts = counts)
           TP10K = norm_TP10K(obj)
           logTP10K = norm_logTP10K(obj)
           SCT = SCTransform(obj, return.only.var.genes = F, verbose = F)
           
           ## norm counts into DT ------------------------------------------------------
           DT_with.meta = function(counts){
             dt.t.counts = t(counts) %>% as.data.table()
             dt.t.counts$CellType = gsub("\\..*", "", colnames(counts))
             dt.t.counts$condition = gsub(".*\\.set", "CD", colnames(counts))
             dt.t.counts$condition = gsub("\\..*", "", dt.t.counts$condition)
             
             return(dt.t.counts)
           }
           
           counts_dt.meta = DT_with.meta(counts); rm(obj)
           TP10K_dt.meta = DT_with.meta(TP10K[["RNA"]]@data); rm(TP10K)
           logTP10K_dt.meta = DT_with.meta(logTP10K[["RNA"]]@data); rm(logTP10K)
           SCT_dt.meta = DT_with.meta(SCT[["SCT"]]@scale.data); rm(SCT)
           
           # ## get gene mean group by cell type ------------------
           # counts_dt.meta.no.cond = counts_dt.meta[,!"condition"]
           # gene.mean = 
           #   counts_dt.meta.no.cond[, lapply(.SD, 
           #                                   function(X){
           #                                     sum(X)/length(X)
           #                                   }), 
           #                          by = c("CellType")] %>% as.data.frame
           # gene.mean = 
           #   reshape2::melt(gene.mean, 
           #                  value.name = "mean") %>% 
           #   dplyr::rename(Gene = variable)
           # 
           # gene.mean$Gene.CellType = paste0(gene.mean$Gene, ".", gene.mean$CellType)
           # gene.mean.altered = gene.mean[gene.mean$Gene.CellType %in% gene.celltype_altered, ]
           # 
           # gene.mean.altered = 
           #   gene.mean.altered %>% group_by(CellType) %>% 
           #   mutate(gene.mean_order = order(mean))
           # 
           # gene.mean.altered = 
           #   gene.mean.altered %>% group_by(CellType) %>% 
           #   mutate(gene.mean_order_quantile = scales::rescale(gene.mean_order))
           # 
           # gene.mean.altered$group = NA
           # gene.mean.altered$group[gene.mean.altered$gene.mean_order_quantile < 0.25] = "Group 4"
           # gene.mean.altered$group[gene.mean.altered$gene.mean_order_quantile >= 0.25 & gene.mean.altered$gene.mean_order_quantile < 0.5] = "Group 3"
           # gene.mean.altered$group[gene.mean.altered$gene.mean_order_quantile >= 0.5 & gene.mean.altered$gene.mean_order_quantile < 0.75] = "Group 2"
           # gene.mean.altered$group[gene.mean.altered$gene.mean_order_quantile >= 0.75 & gene.mean.altered$gene.mean_order_quantile < 1] = "Group 1"
           # 
           # calculate var --------------------------------------------------------
           var_counts = counts_dt.meta[, lapply(.SD, var), by = c("CellType", "condition")] %>% as.data.frame
           var_TP10K = TP10K_dt.meta[, lapply(.SD, var), by = c("CellType", "condition")] %>% as.data.frame
           var_logTP10K = logTP10K_dt.meta[, mclapply(.SD, var, mc.cores = 4), by = c("CellType", "condition")] %>% as.data.frame
           var_SCT = SCT_dt.meta[, mclapply(.SD, var, mc.cores = 4), by = c("CellType", "condition")]%>% as.data.frame
           
           # merge and melt -------------------------------------------------------
           ## MeltVar function
           MeltVar = function(var){
             melt.var = reshape2::melt(var, id = c("CellType", "condition"), value.name = "var")
             return(melt.var)
           }
           ## MeltVar: melt vars
           var_counts.melt = MeltVar(var_counts)
           var_TP10K.melt = MeltVar(var_TP10K)
           var_logTP10K.melt = MeltVar(var_logTP10K)
           var_SCT.melt =  MeltVar(var_SCT)
           var_GTruth = var_copula.merged_altered_reform
           var_GTruth.melt = reshape2::melt(var_GTruth, by = c("CellType", "Gene.CellType"), 
                                            value.name = "var")
           var_GTruth.melt = rename(var_GTruth.melt, variable = Gene.CellType, 
                                    condition = variable)
           var_GTruth.melt$variable = gsub("\\..*", "", var_GTruth.melt$variable)
           ## merge
           ### merge all melt vars
           var_melt.merge = rbind(
             var_counts.melt %>% dplyr::rename(Gene = variable) %>% mutate(NormMethod = "count"), 
             var_TP10K.melt %>% dplyr::rename(Gene = variable) %>% mutate(NormMethod = "TP10K"), 
             var_logTP10K.melt %>% dplyr::rename(Gene = variable) %>% mutate(NormMethod = "logTP10K"),
             var_SCT.melt %>% dplyr::rename(Gene = variable) %>% mutate(NormMethod = "SCT"), 
             var_GTruth.melt %>% dplyr::rename(Gene = variable) %>% mutate(NormMethod = "GT")
           )
           ### merge with gene mean group
           var_melt.merge$Gene.CellType <- paste0(var_melt.merge$Gene, ".", var_melt.merge$CellType)
           # var_melt.merge_by.genegroup = 
           #   merge(var_melt.merge, gene.mean.altered %>% ungroup %>% select(Gene.CellType, mean), 
           #         by = "Gene.CellType", 
           #         all.x = TRUE)
           ### remove NA rows
           # var_melt.merge_by.genegroup = var_melt.merge_by.genegroup[!is.na(var_melt.merge_by.genegroup$group), ]
           
           ### add log var to var_melt.merge
           var_melt.merge$logVar = log10(var_melt.merge$var + 1)
           
           ### add celltype gene
           var_melt.merge$Gene.CellType = paste0(var_melt.merge$Gene, ".", var_melt.merge$CellType)
           
           ## subset only altered genes
           var_melt.merge = var_melt.merge[var_melt.merge$Gene.CellType %in% gene.celltype_altered,  ]
           
           ## get var median --------------------------------------------
           var.median_counts.melt.altered  = 
             setDT(var_melt.merge)[, .(var.median=median(var)), 
                                                .(CellType, condition, NormMethod)]
           ## dcast var mean 
           var.median_counts.melt.altered = dcast(var.median_counts.melt.altered, 
                                                  CellType + NormMethod ~ condition, 
                                                  value.var = 'var.median'
           )
           var.median_counts.melt.altered$diff_var = var.median_counts.melt.altered$CD2 - var.median_counts.melt.altered$CD1
           var.median_counts.melt.altered$FC_var = var.median_counts.melt.altered$CD2/var.median_counts.melt.altered$CD1
           var.median_counts.melt.altered$logFC_var = log2(var.median_counts.melt.altered$FC_var)
           

           ## add batch and nCells information ---------------------------------
           var.median_counts.melt.altered$batch = X
           var.median_counts.melt.altered$ncells_per_type = ncells_per_type
           var.median_counts.melt.altered$ncells_per_type.char = paste0(ncells_per_type, " cells/type")
           
           # writing out -------------------------------------------------------
           dir.create(
             sprintf("XXXX/out/%s-var%.0f_%.1f_%.1f-6type-%scell", 
                     input_datetag, gene_perc, factor_lower, factor_upper, 
                     ncells)
           )
           
           dir.create(
             sprintf("XXXX/out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                     input_datetag, gene_perc, factor_lower, factor_upper, 
                     ncells, X)
           )
           
           outpath = sprintf("XXXX/out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                             input_datetag, gene_perc, factor_lower, factor_upper, 
                             ncells, X)
           
           # print("Saving files")
           write.table(as.data.frame(var.median_counts.melt.altered), 
                       #"%s/MedianVar.txt", 
                       sprintf("%s/Median.Var.txt", outpath), 
                       quote = F, col.names = T, row.names = F, sep = "\t")
           
           # print("finished")           
           print(sprintf("saved output of batch %1.0f of %scells/type, at %s", 
                         X, 
                         ncells, 
                         Sys.time() %>% as.character))
         })
