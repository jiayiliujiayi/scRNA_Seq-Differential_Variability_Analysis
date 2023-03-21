.libPaths(c(.libPaths(), "/home/XXXX/R/rocker-rstudio/4.0/"))

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
library(Matrix, lib.loc = "/home/XXXX/R/rocker-rstudio/4.0/")
dyn.load('/home/XXXX/miniconda3/lib/libxml2.so.2'); dyn.load('/home/XXXX/miniconda3/lib/libglpk.so.40');library(Seurat)
library(magrittr)

datetag = gsub("\\-", "", Sys.Date())

# extract args
args = commandArgs(trailingOnly=T)
folder = args[1] %>% as.character()

minfolder = input_datetag = args[grep('--minfolder=', args)] %>% gsub("--minfolder=", "", .)
maxfolder = input_datetag = args[grep('--maxfolder=', args)] %>% gsub("--maxfolder=", "", .)
input_datetag = args[grep('--input_datetag=', args)] %>% gsub("--input_datetag=", "", .) %>% as.character
gene_perc = args[grep('--gene_perc=', args)] %>% gsub("--gene_perc=", "", .) %>% as.numeric
factor_upper = args[grep('--factor_upper=', args)] %>% gsub("--factor_upper=", "", .) %>% as.numeric
factor_lower = args[grep('--factor_lower=', args)] %>% gsub("--factor_lower=", "", .) %>% as.numeric
ncells_per_type = args[grep('--ncells_per_type=', args)] %>% gsub("--ncells_per_type=", "", .) %>% as.numeric
# HERE Change n Types
ncells = 6 * ncells_per_type


mclapply(minfolder:maxfolder, mc.cores = 8, 
         function(X){
           path = sprintf("/home/XXXX/denSNE/8.0calculate-affinity_batches/out/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                          input_datetag, 
                          gene_perc, 
                          factor_lower, 
                          factor_upper, 
                          ncells, X)
           print(sprintf("Started wilcox test for dataset %1.0f of %scells/type at %s", 
                         X, 
                         ncells, 
                         as.character(Sys.time())))
           
           # import ------------------------------------------------------------------
           files = list.files(path = path, pattern = "method", full.names = F)
           file.names = gsub(".*\\-", "", files) %>% gsub(".txt", "", .)

           affinities_list = lapply(files,  
                                    function(X){
                                      dat = read.delim(paste0(path, "/", X))
                                      colnames(dat) = c("Cell", gsub(".*\\-", "", X) %>% gsub(".txt", "", .))
                                      return(dat)})
           names(affinities_list) = file.names
           
           # metadata -------------------------------------------------------------
           metadata = data.frame(Cell = affinities_list$method1$Cell, 
                                 CellType =  affinities_list$method1$Cell %>% gsub("\\..*", "",.), 
                                 sim.method = affinities_list$method1$Cell %>% gsub("[^\\.]*\\.(.*)", "\\1", .))
           metadata$sim.method = gsub("\\..*", "", metadata$sim.method)

           # merge metadata ----------------------------------------------------------
           meta.affinities_list = 
             lapply(affinities_list, function(X) merge(X, metadata))

           # stat test ---------------------------------------------------------------
           p.w_list = 
             mclapply(meta.affinities_list, mc.cores = 8, 
                      function(X){
               # init p d values data frame
               p_w_values = data.frame(CellType = unique(metadata$CellType) %>% sort, 
                                       affinity = colnames(X)[2], 
                                       W = NA, 
                                       p.val = NA, 
                                       median.diff_2to1 = NA
               )
               
               for (celltype in unique(p_w_values$CellType)){
                 X.sub = X[X$CellType == celltype, ]
                 aff.set1 = X.sub[X.sub$sim.method == "set1", 2]
                 aff.set2 = X.sub[X.sub$sim.method == "set2", 2]
                 
                 if (celltype %in% c("B cells", "Monocytes")) {
                   tryCatch(
                     {
                       out.mww = wilcox.test(x = aff.set1, y = aff.set2, 
                                             alternative = "less", 
                                             paired = F)
                     }, 
                     error = function(e) {}
                   )
                 } else if (celltype %in% c("Neutrophils", "Platelets")){
                   tryCatch(
                     {
                       out.mww = wilcox.test(x = aff.set1, y = aff.set2, 
                                             alternative = "greater", 
                                             paired = F)
                     }, 
                     error = function(e) {}
                   )
                 } else if (celltype %in% c("RBC", "T cells")) {
                   tryCatch(
                     {
                       out.mww = wilcox.test(x = aff.set1, y = aff.set2, 
                                             alternative = "two.sided", 
                                             paired = F)
                     }, 
                     error = function(e) {}
                   )
                   
                 }
                 
                 
                 tryCatch(
                   {
                     p_w_values[p_w_values$CellType == celltype, "W"] = out.mww$statistic
                     p_w_values[p_w_values$CellType == celltype, "p.val"] = out.mww$p.value
                     p_w_values[p_w_values$CellType == celltype, "median.diff_2to1"] = median(aff.set2) - median(aff.set1)
                   }, 
                   error = function(e) {}
                 )
                 
               }
               
               return(p_w_values)
               
             })
           
           # violinplot --------------------------------------------------------------
           plot_list = 
             mclapply(meta.affinities_list, mc.cores = 8,
                      function(X) {
               dat = X
               aff.method = colnames(dat)[2]
               p.vals = p.w_list[[aff.method]]
               p.vals$sim.method = "none"
               colnames(dat)[2] = "affinity"
               p.vals$sig = ifelse(p.vals$p.val < 0.05, "sig", "non.sig")
               p.vals$label = paste0("p.val=", as.character(format(p.vals$p.val, scientific = T, digits = 2)), 
                                     ", ", p.vals$sig)
               ggplot(dat, aes(x = sim.method, y = affinity)) + 
                 geom_violin(aes(color = sim.method)) + 
                 geom_boxplot(outlier.alpha = 0, aes(color = sim.method), width = 0.4) + 
                 labs(x = "sim group", y = paste0("affinities from ", aff.method)) +
                 facet_wrap( ~ CellType, scales = "free", nrow = 3) +
                 scale_color_manual(values = c("gray", "firebrick3")) +
                 ggtitle(aff.method) + 
                 theme_minimal() + 
                 geom_label(data = p.vals, aes(label = label), inherit.aes = FALSE, 
                            x = Inf, y = Inf, hjust=1, vjust=1
                 )
             })
           
           
           # save files --------------------------------------------------------------
           dir.create(
             sprintf("./out.MWW/%s-var%.0f_%.1f_%.1f-6type-%scell", 
                     input_datetag, gene_perc, factor_lower, factor_upper, 
                     ncells)
           )
           
           dir.create(
             sprintf("./out.MWW/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                     input_datetag, gene_perc, factor_lower, factor_upper, 
                     ncells, X)
           )
           
           outpath = sprintf("./out.MWW/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
                             input_datetag, gene_perc, factor_lower, factor_upper, 
                             ncells, X)
           
           saveRDS(plot_list, sprintf("%s/list-plot.Rds", outpath))
           saveRDS(p.w_list, sprintf("%s/list-p_w.val.Rds", outpath))
           saveRDS(meta.affinities_list, sprintf("%s/list-affinities.Rds", outpath))
           
           print(sprintf("saved results for dataset %1.0f of %scells/type, at %s", 
                         X, 
                         ncells, 
                         Sys.time() %>% as.character))
           }

         
         )



# # import ------------------------------------------------------------------
# ## list files
# print("fetching affinity files")
# files = list.files(path = path, pattern = "method", full.names = F)
# file.names = gsub(".*\\-", "", files) %>% gsub(".txt", "", .)
# print(files)
# # if(length(files) == 11){
# #   print("all affinity files are available")
# # }else if (length(files) != 11) {
# #   sprintf("only %d files are available", length(files))
# #   # get missing files
# #   files.num = gsub("method", "", file.names) %>% as.numeric
# #   notaval.files = unique(files.num, 1:11)
# #   print("missing affinities from method(s): "); print(notaval.files)
# # }
# 
# print("importing datasets")
# affinities_list = lapply(files, 
#                          function(X){
#                            dat = read.delim(paste0(path, "/", X))
#                            colnames(dat) = c("Cell", gsub(".*\\-", "", X) %>% gsub(".txt", "", .))
#                            return(dat)})
# names(affinities_list) = file.names
# 
# # metadata -------------------------------------------------------------
# print("Fetching metadata")
# metadata = data.frame(Cell = affinities_list$method1$Cell, 
#                       CellType =  affinities_list$method1$Cell %>% gsub("\\..*", "",.), 
#                       sim.method = affinities_list$method1$Cell %>% gsub("[^\\.]*\\.(.*)", "\\1", .))
# metadata$sim.method = gsub("\\..*", "", metadata$sim.method)
# 
# 
# # merge metadata ----------------------------------------------------------
# meta.affinities_list = 
#   lapply(affinities_list, function(X) merge(X, metadata))
# 
# 
# # stat test ---------------------------------------------------------------
# print("Performing stat test")
# 
# p.w_list = 
# lapply(meta.affinities_list, function(X){
#   # init p d values data frame
#   p_w_values = data.frame(CellType = unique(metadata$CellType) %>% sort, 
#                           affinity = colnames(X)[2], 
#                           W = NA, 
#                           p.val = NA, 
#                           median.diff_2to1 = NA
#                           )
#   
#   for (celltype in unique(p_w_values$CellType)){
#     X.sub = X[X$CellType == celltype, ]
#     aff.set1 = X.sub[X.sub$sim.method == "set1", 2]
#     aff.set2 = X.sub[X.sub$sim.method == "set2", 2]
#     
#     if (celltype %in% c("B cells", "Monocytes")) {
#       tryCatch(
#         {
#           out.mww = wilcox.test(x = aff.set1, y = aff.set2, 
#                                 alternative = "less", 
#                                 paired = F)
#         }, 
#         error = function(e) {}
#       )
#     } else if (celltype %in% c("Neutrophils", "Platelets")){
#       tryCatch(
#         {
#           out.mww = wilcox.test(x = aff.set1, y = aff.set2, 
#                                 alternative = "greater", 
#                                 paired = F)
#         }, 
#         error = function(e) {}
#       )
#     } else if (celltype %in% c("RBC", "T cells")) {
#       tryCatch(
#         {
#           out.mww = wilcox.test(x = aff.set1, y = aff.set2, 
#                                 alternative = "two.sided", 
#                                 paired = F)
#         }, 
#         error = function(e) {}
#       )
#       
#     }
#     
# 
#     tryCatch(
#       {
#         p_w_values[p_w_values$CellType == celltype, "W"] = out.mww$statistic
#         p_w_values[p_w_values$CellType == celltype, "p.val"] = out.mww$p.value
#         p_w_values[p_w_values$CellType == celltype, "median.diff_2to1"] = median(aff.set2) - median(aff.set1)
#       }, 
#       error = function(e) {}
#     )
#     
#   }
#   
#   #p_w_values$direction[p_w_values$CellType %in% c("B cells", "Monocytes")] = p_w_values$median.diff_2to1[p_w_values$CellType %in% c("B cells", "Monocytes")] > 0
#   #p_w_values$direction[p_w_values$CellType %in% c("Neutrophils", "Platelets")] = p_w_values$median.diff_2to1[p_w_values$CellType %in% c("Neutrophils", "Platelets")] < 0
#   
#   
#   return(p_w_values)
#   
# })
# 
# # violinplot --------------------------------------------------------------
# plot_list = 
#   lapply(meta.affinities_list, function(X) {
#     dat = X
#     aff.method = colnames(dat)[2]
#     p.vals = p.w_list[[aff.method]]
#     p.vals$sim.method = "none"
#     colnames(dat)[2] = "affinity"
#     p.vals$sig = ifelse(p.vals$p.val < 0.05, "sig", "non.sig")
#     p.vals$label = paste0("p.val=", as.character(format(p.vals$p.val, scientific = T, digits = 2)), 
#                           ", ", p.vals$sig)
#     ggplot(dat, aes(x = sim.method, y = affinity)) + 
#       geom_violin(aes(color = sim.method)) + 
#       geom_boxplot(outlier.alpha = 0, aes(color = sim.method), width = 0.4) + 
#       labs(x = "sim group", y = paste0("affinities from ", aff.method)) +
#       facet_wrap( ~ CellType, scales = "free", nrow = 3) +
#       scale_color_manual(values = c("gray", "firebrick3")) +
#       ggtitle(aff.method) + 
#       theme_minimal() + 
#       geom_label(data = p.vals, aes(label = label), inherit.aes = FALSE, 
#                  x = Inf, y = Inf, hjust=1, vjust=1
#                 )
#   })
# 
# 
# # save files --------------------------------------------------------------
# print("Saving files")
# 
# dir.create(
#   sprintf("./out.MWW/%s-var%.0f_%.1f_%.1f-6type-%scell", 
#           input_datetag, gene_perc, factor_lower, factor_upper, 
#           ncells)
# )
# 
# dir.create(
#   sprintf("./out.MWW/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
#           input_datetag, gene_perc, factor_lower, factor_upper, 
#           ncells, folder)
# )
# 
# outpath = sprintf("./out.MWW/%s-var%.0f_%.1f_%.1f-6type-%scell/%s", 
#                   input_datetag, gene_perc, factor_lower, factor_upper, 
#                   ncells, folder)
# 
# saveRDS(plot_list, sprintf("%s/list-plot.Rds", outpath))
# saveRDS(p.w_list, sprintf("%s/list-p_w.val.Rds", outpath))
# saveRDS(meta.affinities_list, sprintf("%s/list-affinities.Rds", outpath))
# 
# # finishing ---------------------------------------------------------------
# print("Finished")
# 
