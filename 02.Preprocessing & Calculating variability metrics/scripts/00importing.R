importData = function(path){
  # import raw data ---------------------------------------------------------
  counts.set1 = readRDS(paste0(path, "/set1.rds"))
  counts.set2 = readRDS(paste0(path, "/set2.rds"))
  testcounts = readRDS(paste0("~/denSNE/6.0simulation/", "20211014", "-test_counts_multi_type.blood.rds"))
  rownames(counts.set1) = rownames(testcounts)
  rownames(counts.set2) = rownames(testcounts)
  
  # init Seurat obj ---------------------------------------------------------
  obj.set1 = CreateSeuratObject(counts = counts.set1)
  obj.set1$Major.cell.type = colnames(counts.set1) %>% gsub("\\..*", "",.)
  obj.set1$sim.method = colnames(counts.set1) %>% gsub("^.*\\.*", "",.)
  
  obj.set2 = CreateSeuratObject(counts = counts.set2)
  obj.set2$Major.cell.type = colnames(counts.set2) %>% gsub("\\..*", "",.)
  obj.set2$sim.method = colnames(counts.set2) %>% gsub("^.*\\.", "",.)
  
  
  # merge objs --------------------------------------------------------------
  obj = merge(obj.set1, y = obj.set2)
  obj$sim.group = ifelse(obj$sim.method == "set1", "set1", "set2")
  obj$sim.method = factor(obj$sim.method, levels = c("set1", "set2"))
  obj$Cell = colnames(obj)
  
  # return metadata, counts -------------------------------------------------
  metadata = obj@meta.data
  counts.raw = obj[["RNA"]]@counts
  
  out.list = list(obj = obj, 
                  metadata = metadata, 
                  counts.raw = counts.raw)
  
  return(out.list)
}


