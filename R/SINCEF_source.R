#'@title process function
#'
#'@description An unsupervised clustering approach for Single-cell methylome clustering by reconstructing dissimilarity matrix between pairs of cells with embedding fusion.
#'
#'@param dataframe
#'
#'@return dataframe
#'
#'@examples df_test <- mat_F(df)
#'
#'@export
mat_F <- function(df){
  chr_P <- paste(df[,1], df[,2], sep = "_")
  df <- data.frame(chr_P,df)
  df <- df[,c(1,4)]
  return(df)
  rm(chr_P)
}

#'@title Process function
#'
#'@description Get the corresponding dissimilarity matrix using single-cell methylation data
#'
#'@param list
#'
#'@return Pearson dissimilarity object
Get_Pearson <- function(j,data){

  cor_cp <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    pv <- cor(dm[,2],dm[,3])
    return(pv)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate Pearson correlation coefficient

  res_cor <- lapply(1:length(data),cor_cp,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@title Process function
#'
#'@description Get the corresponding dissimilarity matrix using single-cell methylation data
#'
#'@param list
#'
#'@return Cosine dissimilarity object
Get_Cosine <- function(j,data){

  cor_cos <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    cos <- sum(dm[,2]*dm[,3])/sqrt((sum(dm[,2]^2)*sum(dm[,3]^2)))
    return(cos)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate cosine correlation coefficient

  res_cor <- lapply(1:length(data),cor_cos,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@title Process function
#'
#'@description Get the corresponding dissimilarity matrix using single-cell methylation data
#'
#'@param list
#'
#'@return Dual dissimilarity object
Get_Dual <- function(j,data){

  cor_dua <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    dm$cot <- ifelse(dm[,2]==dm[,3],1,0)
    dua <- sum(dm$cot)/length(dm$cot)
    return(dua)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate cosine correlation coefficient

  res_cor <- lapply(1:length(data),cor_dua,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@title Dissimilarity object
#'
#'@description Get the corresponding dissimilarity matrix using single-cell methylation data
#'
#'@param k_cpu number of CPU cores for parallel computing
#'
#'@param data R list which is composed of dataframes of cell methylation information, including "Chromosome number", "locus position", and "methylation_level(0 or 1)"
#'
#'@param method the distance measure to be used. This must be one of "Cosine" ,"Dual" or "Pearson"
#'
#'@return the corresponding dissimilarity matrix
#'
#'@examples
#'data_cell <- load('data_blood.RData')
#'view(data_cell)
#'dism_pearson <- Output_DISM(k_cpu=8,data_cell,method='Pearson')
#'
#'@export
Output_DISM <- function(k_cpu,data,method){
  if(method=='Cosine'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Cosine,data)
    stopCluster(cl)
    res_mat <- do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  if(method=='Dual'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Dual,data)
    stopCluster(cl)
    res_mat<-do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  if(method=='Pearson'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Pearson,data)
    stopCluster(cl)
    res_mat<-do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  else{
    print('Something wrong with your settings,please check...')
  }

}

#'@title Normalization
#'
#'@description Normalize a vector or matrix
#'
#'@param x vectors or matrix
#'
#'@param method the normalization method to be used. This must be one of "Norm" ,or "scale"
#'
#'@return vectors or matrix
#'
#'@examples x_nm <- nm(x,method='Norm')
#'
#'@export
nm<-function(x){
    num <- ncol(x)
    for (i in 1:num) {
      x[,i] <- (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))
    }
    return(x)
}

#'@title Dissimilarity matrix reconstruction
#'
#'@description The input distance matrices are embedded and fused, and the normalized reconstructed distance matrix is output
#'
#'@param d1 one of the three dissimilarity matrices, including Pearson_dism, Cosine_dism and Dual_dism
#'
#'@param d2 one of the three dissimilarity matrices, including Pearson_dism, Cosine_dism and Dual_dism
#'
#'@param d3 one of the three dissimilarity matrices, including Pearson_dism, Cosine_dism and Dual_dism
#'
#'@param nei_k the n_neighbors to be used in UMAP transformation. The recommended value range is 5~30.
#'
#'@param is_scale whether to use normalization for reconstructing dissimilarity matrix.
#'
#'@return List which is composed of normalized reconstructed dissimilarity matrix,and the other three expanded dissimilarity matrices
#'
#'@examples
#'dism_pearson <- Output_DISM(k_cpu=8,data_cell,method='Pearson')
#'dism_cosine  <- Output_DISM(k_cpu=8,data_cell,method='Cosine')
#'dism_dual    <- Output_DISM(k_cpu=8,data_cell,method='Dual')
#'rsm          <- Get_resm(dism_cosine,dism_dual,dism_pearson,6,is_scale="TRUE")
#'
#'@export
Get_resm <- function(d1,d2,d3,nei_k,is_scale){
  if(is_scale=='TRUE'){
    f1  <- umap(d1,input='dist',n_neighbors=nei_k)
    f2  <- umap(d2,input='dist',n_neighbors=nei_k)
    f3  <- umap(d3,input='dist',n_neighbors=nei_k)
    em1 <- nm(f1$layout)
    em2 <- nm(f2$layout)
    em3 <- nm(f3$layout)
    em  <- as.matrix(dist(cbind(em1,em2,em3)))
    return(list(em,em1,em2,em3))
  }
  if(is_scale=='FALSE'){
    f1  <- umap(d1,input='dist',n_neighbors=nei_k)
    f2  <- umap(d2,input='dist',n_neighbors=nei_k)
    f3  <- umap(d3,input='dist',n_neighbors=nei_k)
    em1 <- f1$layout
    em2 <- f2$layout
    em3 <- f3$layout
    em  <- as.matrix(dist(cbind(em1,em2,em3)))
    return(list(em,em1,em2,em3))
  }
  else{
    print('Check your settings...')
  }

}

#'@title Hierarchical clustering
#'
#'@description Input the reconstructed dissimilarity matrix and output the hierarchical clustering result to help judge the number of potential cell clusters.
#'
#'@param rsm_list the reconstruction dissimilarity matrix list
#'
#'@examples Get_hclust(rsm)
#'
#'@return  hierarchical clustering result
#'
#'@export
Get_hclust <- function(rsm_list){
  return(pheatmap(rsm_list[[1]],cutree_rows = NA,cutree_cols=NA,
                  cellwidth = 5, cellheight = 5))
}

#'@title Assess clustering performace
#'
#'@description Calculate ARI to evaluate clustering performance in combination with known cell types
#'
#'@param rsm_list the reconstruction dissimilarity matrix list
#'
#'@param clu the potential cluster number which could be obtained from hierarchical clustering heatmap of Get_hclust
#'
#'@param ref_label the available or known cell type label
#'
#'@examples
#'true_label <- load('syth_label.RData')
#'k_clu <- 6
#'ARI_test <- Comp_ARI(rsm,k_clu,true_label)
#
#'@return the clustering performance
#'
#'@export
Comp_ARI <- function(rsm_list,clu,ref_label){
  hclu   <- hclust(dist(cbind(rsm_list[[2]],rsm_list[[3]],rsm_list[[4]])))
  pre    <- cutree(hclu,k=clu)
  return(ARI(pre,ref_label))
}
