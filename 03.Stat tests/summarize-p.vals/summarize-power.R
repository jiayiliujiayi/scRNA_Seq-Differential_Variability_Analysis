library(GDAtools)
library(rdist)
library(dplyr)     
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggplot2)
library(Matrix, lib.loc = "/home/XXXX/R/rocker-rstudio/4.0/")
dyn.load('/home/XXXX/miniconda3/lib/libxml2.so.2')
dyn.load("/home/XXXX/miniconda3/lib/libXt.so.6")
dyn.load('/home/XXXX/miniconda3/lib/libglpk.so.40');library(Seurat)
library(magrittr)
library(corrplot)
library(RColorBrewer)

datetag = gsub("-", "", Sys.Date())

# calc power for MWW
paths = list.dirs(path="/home/XXXX/denSNE/9.0stat_batches/out.MWW", recursive = F) %>% 
  .[grep("20220309", .)]
paths = paste0(paths, "/")

calc.power = 
  function(path){
    files = list.files(path, pattern = "list-p_w.val.Rds", recursive = T)
    # import as list
    p_w.list = 
      lapply(files, 
             function(X){
               file.path = paste0(path, X)
               dat = readRDS(file.path)
               return(dat)
             }
      )
    
    names(p_w.list) = gsub("/list-p_w.val.Rds", "", files) %>% paste0("batch", .)
    
    p_w.list.clean = p_w.list
    for (i in 1:length(p_w.list.clean)) {
      each.pair_list = p_w.list.clean[[i]]
      for (each.method in names(each.pair_list)) {
        each.pair_list[[each.method]] = mutate(each.pair_list[[each.method]], method = each.method)
      }
      p_w.list.clean[[i]] = each.pair_list
    }
    
    # transform list of lists into list of data frames
    p_w.list.clean = lapply(p_w.list.clean, function(X) Reduce(rbind, X))
    
    # mutate batch name
    p_w.list.clean = lapply(names(p_w.list.clean), 
                            function(X)
                              p_w.list.clean[[X]] %>% mutate(batch = X))
    # transform list of dataframes into dataframe
    p_w = Reduce(rbind, p_w.list.clean)
    
    p_w = 
      p_w %>%
      group_by(CellType, method) %>% 
      mutate(p.adj = p.adjust(p.val, method='fdr'))
    
    p_w$sig = ifelse(p_w$p.adj < 0.05, "sig", "non.sig")
    
    # calculate power
    power = p_w %>% select(CellType, method, sig) %>% group_by(CellType, method) %>% summarise(power = length(which(sig == "sig"))/length(sig))
    power = dcast(power, CellType ~ method)
    rownames(power) = power$CellType
    power = as.matrix(power[, -1])
    
    return(power)
  }

list.power = lapply(paths, calc.power)

names(list.power) = gsub("/home/XXXX/denSNE/9.0stat_batches/out.MWW/20220309-var50_1.5_3.0-6type-", "", paths)
names(list.power) = gsub("\\/", "", names(list.power))

saveRDS(list.power, 
        sprintf("/home/XXXX/denSNE/9.0stat_batches/summarize-p.vals/%s-list-power-MWW.Rds", datetag))


