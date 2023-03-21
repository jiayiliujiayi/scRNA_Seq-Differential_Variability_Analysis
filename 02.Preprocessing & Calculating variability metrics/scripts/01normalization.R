# TP10K -------------------------------------------------------------------
norm_TP10K <- function(obj){
  # --------- depcrecated part, use Seurat::Normalize(method = "RC") instead------------
  counts.raw = as.matrix(obj[["RNA"]]@counts)

  # calculate lib size
  libsize = colSums(counts.raw)
  # calculate
  tp10k = sweep(counts.raw, 2, libsize, FUN = '/')*1e4
  tp10k = as.matrix(tp10k)

  #tp10k = as(tp10k, "dgCMatrix")

  obj.new <- SeuratObject::SetAssayData(object = obj, assay = "RNA", slot = "data",
                                        new.data = tp10k)
  # --------- depcrecated part, use Seurat::Normalize(method = "RC") instead------------
  #obj.new = NormalizeData(obj, normalization.method = "RC")
  
  return(obj.new)
}


# logTP10K ----------------------------------------------------------------
norm_logTP10K <- function(obj){
  counts.raw = as.matrix(obj[["RNA"]]@counts)
  
  # calculate lib size
  libsize = colSums(counts.raw)
  # calculate
  tp10k = sweep(counts.raw, 2, libsize, FUN = '/')*1e4
  logtp10k = log(tp10k + 1)
  
  obj.new <- SetAssayData(obj, assay = "RNA", slot = "data", 
                      new.data = as.matrix(logtp10k))
  
  return(obj.new)
}

# scTransform -------------------------------------------------------------
norm_SCT <- function(obj){
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  return(obj)
}
