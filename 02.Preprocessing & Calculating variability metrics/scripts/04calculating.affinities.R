source("~/pkg/densne/densne.R")
calc_affinity = function(obj, method = NULL, cell.type = NULL){
  if (method == "local.radius" && is.null(cell.type)){
    coord_pca = Embeddings(obj, reduction = "pca")
    coord_pca = coord_pca[, 1:min(200, ncol(coord_pca))]
    out_densne = run_densne(coord_pca, no_dims=2, perplexity=50,
                            theta=0.5,
                            randseed=-1,
                            verbose=FALSE, use_pca=FALSE,
                            max_iter=1000, dens_frac=0.3, # as in the paper
                            final_dens=TRUE)
    local.radius = out_densne[[2]]

    out = data.frame(Cell = colnames(obj),
                     aff = local.radius)
  } else if (method == 'embedding.to.LR_ro' && is.null(cell.type)){
    embeddings = Embeddings(obj, reduction = "denSNE")
    out_densne = run_densne(embeddings, no_dims=2, perplexity=50,
                            theta=0.5,
                            randseed=-1,
                            verbose=FALSE, use_pca=FALSE,
                            max_iter=1000, dens_frac=0.3, # as in the paper
                            final_dens=TRUE)
    local.radius = out_densne[[2]]
    
    out = data.frame(Cell = colnames(obj),
                     aff = local.radius)
  } else if (method == 'local.radius_re' && is.null(cell.type)){
    coord_pca = Embeddings(obj, reduction = "pca")
    coord_pca = coord_pca[, 1:min(200, ncol(coord_pca))]
    out_densne = run_densne(coord_pca, no_dims=2, perplexity=50,
                            theta=0.5,
                            randseed=-1,
                            verbose=FALSE, use_pca=FALSE,
                            max_iter=1000, dens_frac=0.3, # as in the paper
                            final_dens=TRUE)
    local.radius = out_densne[[3]]
    
    out = data.frame(Cell = colnames(obj),
                     aff = local.radius)
  } else if (method == "dist.medoid" && !is.null(cell.type)) {
    library(GDAtools)
    library(rdist)
    # get embeddings and cell types
    coord_pca = Embeddings(obj, reduction = "pca")
    clust = obj@meta.data[, c("sim.method", cell.type)]
    clust$Cell = colnames(obj)
    
    # calculate distances
    dist_pca = stats::dist(coord_pca, method = "euclidean")
    
    # finding medoids
    ## get cluster vector
    clust$cell_type.sim_method = paste0(clust$Major.cell.type, ".", clust$sim.method)
    clust$cell_type.sim_method_cluster = as.numeric(as.factor(clust$cell_type.sim_method))
    ## identifying medoids
    medoid = medoids(dist_pca, clust$cell_type.sim_method_cluster)
    names(medoid) = clust$cell_type.sim_method %>% unique %>% sort
    ## get medoids coordinates
    coord_medoid_pcs = coord_pca[medoid, ]
    rownames(coord_medoid_pcs) = names(medoid)
    
    # calculating each distances
    medoid_dist = rep(NA, nrow(coord_pca))
    names(medoid_dist) = clust$Cell
    for (cellid in clust$Cell) {
      # get cluster
      cluster = clust$cell_type.sim_method[clust$Cell == cellid]
      
      ## get PCs of cell and medoid
      cell_pcs = coord_pca[which(clust$Cell == cellid), ]
      medoid_pcs = coord_medoid_pcs[cluster, ]
      
      ## get distances and write out
      medoid_dist[cellid] = cdist(cell_pcs, medoid_pcs)[1, 1]
      
      # out data frame
      out = data.frame(Cell = names(medoid_dist), aff = medoid_dist, 
                       stringsAsFactors = F, check.names = F)
    }
  } else if (method == "2d.denSNE_to_DM" && !is.null(cell.type)) {
    library(GDAtools)
    library(rdist)
    # get embeddings and cell types
    coord_pca = Embeddings(obj, reduction = "denSNE")
    clust = obj@meta.data[, c("sim.method", cell.type)]
    clust$Cell = colnames(obj)
    
    # calculate distances
    dist_pca = stats::dist(coord_pca, method = "euclidean")
    
    # finding medoids
    ## get cluster vector
    clust$cell_type.sim_method = paste0(clust$Major.cell.type, ".", clust$sim.method)
    clust$cell_type.sim_method_cluster = as.numeric(as.factor(clust$cell_type.sim_method))
    ## identifying medoids
    medoid = medoids(dist_pca, clust$cell_type.sim_method_cluster)
    names(medoid) = clust$cell_type.sim_method %>% unique %>% sort
    ## get medoids coordinates
    coord_medoid_pcs = coord_pca[medoid, ]
    rownames(coord_medoid_pcs) = names(medoid)
    
    # calculating each distances
    medoid_dist = rep(NA, nrow(coord_pca))
    names(medoid_dist) = clust$Cell
    for (cellid in clust$Cell) {
      # get cluster
      cluster = clust$cell_type.sim_method[clust$Cell == cellid]
      
      ## get PCs of cell and medoid
      cell_pcs = coord_pca[which(clust$Cell == cellid), ]
      medoid_pcs = coord_medoid_pcs[cluster, ]
      
      ## get distances and write out
      medoid_dist[cellid] = cdist(cell_pcs, medoid_pcs)[1, 1]
      
      # out data frame
      out = data.frame(Cell = names(medoid_dist), aff = medoid_dist, 
                       stringsAsFactors = F, check.names = F)
    }
  }
  
  return(out)
}
