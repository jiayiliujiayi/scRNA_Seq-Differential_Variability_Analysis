scaleData <- function(obj, cell.type=NULL, ncore=10){
  # function to center and scale data
  func.center_scale = function(X){
    X = (X - mean(X)) / sd(X) ^ as.logical(sd(X))
    return(X)
  }
  
  if(is.null(cell.type)){
    #print("Centering and scaling")
    norm.data = as.matrix(obj[["RNA"]]@data)
    plan("multisession", 
         workers = min(ncore, detectCores())); scale.norm.data = apply(norm.data, 1, 
                                                               func.center_scale)
    
    scale.norm.data = t(scale.norm.data)

    obj <- SetAssayData(obj, assay = "RNA", slot = "scale.data", 
                        new.data = scale.norm.data)
        
  }else if(!is.null(cell.type)){
    #print("Centering and scaling")
    
    # scale data by cell type
    norm.data = as.matrix(obj[["RNA"]]@data)
    norm.data.celltype = data.frame(CellType = obj@meta.data[, cell.type], t(norm.data), check.names = F)
    norm.data.celltype = setDT(norm.data.celltype)
    plan("multisession", 
         workers = min(ncore, detectCores())); norm.data.celltype = norm.data.celltype[, lapply(.SD, 
                                                                                             func.center_scale),  
                                                                                    by = "CellType"]
    norm.data.celltype = data.frame(norm.data.celltype, check.names = F)
    
    norm.data.celltype.matrix = norm.data.celltype[, !colnames(norm.data.celltype) %in% c("CellType")]
    norm.data.celltype.matrix = as.matrix(norm.data.celltype.matrix)
    norm.data.celltype.matrix = t(norm.data.celltype.matrix)
    colnames(norm.data.celltype.matrix) = colnames(norm.data)
    
    obj <- SetAssayData(obj, assay = "RNA", slot = "scale.data", 
                        new.data = norm.data.celltype.matrix)
    
  } else {
    print("please specify cell types. ")
  }
  
  return(obj)
}

