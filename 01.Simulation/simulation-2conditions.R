.libPaths(c("/home/XXXX/R/rocker-rstudio/4.0/"))

library(scDesign2)
dyn.load('/home/XXXX/miniconda3/pkgs/gsl-2.4-h14c3975_4/lib/libgsl.so.23') # essential for loading gsl and copula package
library(copula)    # corKendall
library(Rtsne)
library(plyr)      # mapvalues
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggplot2)
dyn.load('/home/XXXX/miniconda3/lib/libxml2.so.2'); dyn.load('/home/XXXX/miniconda3/lib/libglpk.so.40');library(Seurat)
library(magrittr)

source("~/pkg/densne/densne.R")

# extract args
args = commandArgs(trailingOnly=T)
seed = args[1]

gene_perc = args[grep('--gene_perc=', args)] %>% gsub("--gene_perc=", "", .) %>% as.numeric
factor_upper = args[grep('--factor_upper=', args)] %>% gsub("--factor_upper=", "", .) %>% as.numeric
factor_lower = args[grep('--factor_lower=', args)] %>% gsub("--factor_lower=", "", .) %>% as.numeric

ncells_per_type = args[grep('--ncells_per_type=', args)] %>% gsub("--ncells_per_type=", "", .) %>% as.numeric

# datetag
datetag = gsub("-", "", Sys.Date())
# load data --------------------------------------------------------------------
obj.raw <- readRDS("/home/XXXX/denSNE/0.0processed/obj.immune.RDS")
Idents(obj.raw) = "Tissue"; obj = subset(obj.raw, idents = "blood"); rm(obj.raw)
Idents(obj) = "cell.type"; obj = subset(obj, 
                                        idents = c("B cells", "Monocytes", "Neutrophils", "Platelets", "RBC", "T cells"))
nCell_per.type = table(obj$Major.cell.type)
index = vector(mode = "list", length = length(nCell_per.type)) %>% `names<-`(names(nCell_per.type))
for (i in 1:length(index)) {
  celltype = names(index)[i]
  if (nCell_per.type[celltype] <= 500) {
    index[[celltype]] = which(obj$Major.cell.type == celltype)
  } else{
    set.seed(41); index[[celltype]] = which(obj$Major.cell.type == celltype) %>% .[sample(1:length(.), 500)]
  } #else {}
}
index = unlist(index) %>% unname
cellid = colnames(obj)[index]; cell_type = obj$cell.type[index] #cell_type %>% table
data_mat.raw <- as.matrix(obj[["RNA"]]@counts)
data_mat.raw <- data_mat.raw[, cellid]
rm(obj)

nGenes_per.cell = read.table("/home/XXXX/denSNE/0.0processed/nGenes_per.cell.txt")
nGenes_per.cell = nGenes_per.cell$V2 %>% `names<-`(nGenes_per.cell$V1)
avr_expr <- read.table("/home/XXXX/denSNE/0.0processed/avr_expr.txt")
avr_expr = avr_expr$V2 %>% `names<-`(avr_expr$V1)

# transform to counts ----------------------------------------------------------
data_mat = data_mat.raw

rownames(data_mat) = names(avr_expr)

gene_zeroexpr = names(which(avr_expr == 0))
data_mat = data_mat[! rownames(data_mat) %in% gene_zeroexpr, ]

nGenes_per.cell = nGenes_per.cell[colnames(data_mat)]
for (gene in nrow(data_mat)) {
  for (cell in ncol(data_mat)) {
    data_mat[gene, cell] <-
      data_mat[gene, cell] * nGenes_per.cell[cell] / avr_expr[gene]
  }
}

data_mat = round(data_mat)

# remove spike-in --------------------------------------------------------------
nonspikes <- which(!grepl("ercc", rownames(data_mat), ignore.case = TRUE))
print(paste("number of spike-ins:", nrow(data_mat)-length(nonspikes)))
data_mat <- data_mat[nonspikes, ,drop = FALSE]

# change colnames into cell type -----------------------------------------------
colnames(data_mat) = cell_type

# split data into train and testing datasets------------------------------------
unique_cell_type <- names(table(colnames(data_mat)))
set.seed(1)
train_idx <- unlist(sapply(unique_cell_type, function(x){
  cell_type_idx <- which(colnames(data_mat) == x)
  n_cell_total <- length(cell_type_idx)
  sample(cell_type_idx, floor(n_cell_total/2))
}))
traincount <- data_mat[, train_idx]
testcount <- data_mat[, -train_idx]

# for reproducibility ----------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")

# select cell types to fit model -----------------------------------------------
nCell_per.type = table(colnames(traincount))
cell_type_sel <- unique(colnames(traincount))
n_cell_new <- ncells_per_type * length(cell_type_sel)
cell_type_prop <- rep(ncells_per_type, 6) %>% `names<-`(cell_type_sel)


# import copula results ---------------------------------------------------
copula_result = readRDS("/home/XXXX/denSNE/6.0simulation/20211014-copula_result_multi_type.blood.default.rds")

# change gene cv: increase variance of percent% of the genes ------------------------
## define a function to change the variance of  genes that are NB fit-----------
change_var <- function(gene_out.list, which_param, 
                       percent.gene, method, 
                       unif.lower, unif.upper, 
                       random.seed){
  # transform X into data frame for easier subsetting
  X = gene_out.list[[which_param]]
  if (!is.null(X)) {
    X = as.data.frame(X)
    colnames(X) = c("percent_zero", "phi", "mu")
    
    # get gene order
    ordered_gene = rownames(X)
    
    # subset genes with NB distribution fitting
    X.nb = X[!is.infinite(X$phi), ]
    
    # randomly select percent% of the genes
    set.seed(random.seed); geneindex = sample(1:nrow(X.nb), 0.01 * percent.gene * nrow(X.nb))
    
    # adding an var column: var = mu + mu^2/phi
    X.nb$var = X.nb$mu + X.nb$mu^2/X.nb$phi
    
    # phi = mu^2 / (var - mu)
    factor = runif(length(geneindex), unif.lower, unif.upper)
    
    # delta
    delta = X.nb$var - X.nb$mu
    
    if (method == "increase") {
      # increase var of the selected genes
      X.nb[geneindex, "var.alt"] = X.nb[geneindex, "var"] * factor
      #X.nb[geneindex, "var"] + 0.999999999 * delta[geneindex]
    } else if (method == "decrease") {
      # decrease var of the selected genes
      var.by.factor = X.nb[geneindex, "var"] / factor
      var.by.mu = 0.01*X.nb[geneindex, "var"] + 0.99 * X.nb[geneindex, "mu"]
      X.nb[geneindex, "var.alt"] = 
        ifelse(
          var.by.factor >  
            var.by.mu, 
          var.by.factor, 
          var.by.mu
        )
    }
    
    # calculate phi.alt
    X.nb[geneindex, "phi"] = X.nb[geneindex, "mu"]^2/(X.nb[geneindex, "var.alt"] - X.nb[geneindex, "mu"])
    
    # replace the NB part in the X
    X[!is.infinite(X$phi),] = X.nb[, -c(4, 5)]
    
    # transform X back
    colnames(X) = NULL
    X = as.matrix(X)
    
    # reorder X
    X = X[ordered_gene, ]
    
    # output gene_out.list
    gene_out.list[[which_param]] = X
  }  else {}
  # return X
  return(gene_out.list)
  
}

inc.id = cell_type_sel[c(1:2)]
dec.id = cell_type_sel[c(3:4)]
rem.id = cell_type_sel[c(5:6)]

# simulate two sets of data
## altering copula results
## set 1: increased variance of decreased cell types, other cell types remain
copula_result.set1 = 
  vector(mode = "list", length = length(cell_type_sel)) %>% `names<-`(cell_type_sel)
### increased variance of decreased cell types
copula_result.set1[dec.id] = 
  lapply(copula_result[dec.id], function(LIST) change_var(LIST, which_param = "marginal_param1", 
                                                          gene_perc, "increase", 
                                                          unif.lower = factor_lower, unif.upper = factor_upper,
                                                          random.seed = seed)
  ) %>% 
  lapply(., function(LIST) change_var(LIST, which_param = "marginal_param2", 
                                      gene_perc, "increase", 
                                      unif.lower = factor_lower, unif.upper = factor_upper, 
                                      random.seed = seed)
  )
### other cell types remain
copula_result.set1[inc.id] = copula_result[inc.id]
copula_result.set1[rem.id] = copula_result[rem.id]
### names to set1
names(copula_result.set1) = paste0(names(copula_result.set1), ".set1")


## set 2: increased variance of increased cell types, other cell types remain
### increased variance of increased cell types
copula_result.set2 = 
  vector(mode = "list", length = length(cell_type_sel)) %>% `names<-`(cell_type_sel)
### increased variance of decreased cell types
copula_result.set2[inc.id] = 
  lapply(copula_result[inc.id], function(LIST) change_var(LIST, which_param = "marginal_param1", 
                                                          gene_perc, "increase", 
                                                          unif.lower = factor_lower, unif.upper = factor_upper,
                                                          random.seed = seed)
  ) %>% 
  lapply(., function(LIST) change_var(LIST, which_param = "marginal_param2", 
                                      gene_perc, "increase", 
                                      unif.lower = factor_lower, unif.upper = factor_upper, 
                                      random.seed = seed)
  )
### other cell types remain
copula_result.set2[dec.id] = copula_result[dec.id]
copula_result.set2[rem.id] = copula_result[rem.id]
### names to set1
names(copula_result.set2) = paste0(names(copula_result.set2), ".set2")


# simulating --------------------------------------------------------------
### gene_perc% [factor_lower, factor_upper]
print(paste0("start simulating data at ", Sys.time()))
sim_counts.set1 <- simulate_count_scDesign2(copula_result.set1, 
                                            n_cell_new, sim_method = 'copula',
                                            cell_type_prop = cell_type_prop %>% `names<-`(names(copula_result.set1)))
sim_counts.set2 <- simulate_count_scDesign2(copula_result.set2, 
                                            n_cell_new, sim_method = 'copula',
                                            cell_type_prop = cell_type_prop %>% `names<-`(names(copula_result.set2)))
print(paste0("finish simulating data at ", Sys.time()))

rownames(sim_counts.set1) = rownames(testcount)
rownames(sim_counts.set2) = rownames(testcount)
# save results ------------------------------------------------------------
# simulation --------------------------------------------------------------
dir.create(
  sprintf("/scratch/XXXX/denSNE/6.1simulation_batches/out/%s-var%.0f_%.1f_%.1f-%stype-%scell", datetag, gene_perc, factor_lower, factor_upper, length(cell_type_sel), n_cell_new)
)

setwd(
  sprintf("/scratch/XXXX/denSNE/6.1simulation_batches/out/%s-var%.0f_%.1f_%.1f-%stype-%scell", datetag, gene_perc, factor_lower, factor_upper, length(cell_type_sel), n_cell_new)
)

dir.create(
  as.character(seed)
)

setwd(
  as.character(seed)
)

saveRDS(sim_counts.set1, 
        file = "set1.rds")
saveRDS(sim_counts.set2, 
        file = "set2.rds")

saveRDS(
  cbind(sim_counts.set1, sim_counts.set2), 
  file = "combined.rds"
)

## copula result
saveRDS(copula_result.set1, 
        file = paste0(datetag, "-","copula_result.set1.rds"))
saveRDS(copula_result.set2, 
        file = paste0(datetag, "-","copula_result.set2.rds"))
