# SINCEF
Single-cell DNA methylome clustering by reconstructing dissimilarity matrix with embedding fusion. SINCEF is provided as an R package which could be installed by the corresponding SINCEF_0.1.0.tar.gz file.

# Tutorial

## Step 1 

   Download SINCEF_0.1.0.tar.gz from the repository and implement local installation via RStudio.
   
   SINCEF relies on the installation of the following R packages for full functionality ："aricode",
                                                                                          "parallel",
                                                                                          "pheatmap",
                                                                                          "umap"
   
   Before using SINCEF, make sure that the above R packages have been installed.
   
## Step 2

   Input your single-cell methylation dataset in the specified format. The recommended format is an R list, where each element is a dataframe containing methylation information for a single cell, organized as follows:
   
    /Chrom/  /location/   /methylation_state/
   
    chr1      131009             1
    
    chr1      131031             0
    
    ...        ...              ...
    
    chrY     56770241            1
    
   More information on data format can be obtained from the simulated dataset "data_bd.RData". 
    
   Download link: https://github.com/TQBio/SINCEF/blob/master/R/data_bd.RData

## Step 3

   Computation of dissimilarity matrices for target single cell methylation data.
    
   Example:
    
    load('data_bd.RData') ###Load simulated data
        
    dism_pearson <- Output_DISM(k_cpu=8,data_BLOOD,method='Pearson')
    
    dism_cosine  <- Output_DISM(k_cpu=8,data_BLOOD,method='Cosine')
    
    dism_dual    <- Output_DISM(k_cpu=8,data_BLOOD,method='Dual')
    
## Step 4

   Using SINCEF to obtain the reconstructed dissimilarity matrix and its corresponding hierarchical clustering object.
   
    test_obj <- Get_resm(dism_cosine, dism_dual, dism_pearson, nei_k = 16, is_scale="TRUE")
    
## Step 5

   Draw a heatmap of hierarchical clustering to identify potential cell subsets.
   
    rcdm <- test_obj[[1]]
    
    Get_hclust(rcdm)

## Step 6

   Load cell reference labels to evaluate hierarchical clustering performance.
   
    load('ref_label.RData')
    
    hc_obj <- test_obj[[2]]
    
    Comp_ARI(hc_obj, clu = 6, cell_bd)
    
# Note

  1. When calculating the dissimilarity matrix, please set it reasonably according to the number of available cores of your PC or server, otherwise it is easy to cause memory overflow.
  
  2. Any problems or bugs encountered during the use of SINCEF, please contact tqglowing@std.uestc.edu.cn.
