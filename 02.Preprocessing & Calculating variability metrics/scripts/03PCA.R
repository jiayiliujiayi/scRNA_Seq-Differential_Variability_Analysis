do_PCA <- function(obj, PCA.input=NULL){
    if (is.null(PCA.input)) {
    print("please specify a PCA method")
      
  }else if (PCA.input == "all.genes") {
    DefaultAssay(obj) = "RNA"
    all.genes = rownames(obj)
    obj = RunPCA(obj, features = all.genes, npcs = min(200, ncol(obj) - 1), verbose = F)
    
  }else if (PCA.input == "variable.genes"){
    DefaultAssay(obj) = "RNA"
    obj = FindVariableFeatures(obj, nfeatures = 2000, verbose = F)
    obj = RunPCA(obj, features = VariableFeatures(obj), npcs = min(200, ncol(obj) - 1), verbose = F)
    
  }else if (PCA.input == "SCT.genes") {
    DefaultAssay(obj) = "SCT"
    obj = RunPCA(obj, npcs = min(200, ncol(obj) - 1), verbose = FALSE)
    
  }
  
  return(obj)
}


