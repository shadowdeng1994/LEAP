#' Title
#'
#' @param SeuratObject A Seurat object.
#' @param nPC Number of PC.
#' @param Verbose Verbose or not.
#'
#' @returns A processed Seurat object
#' @export
#'
#' @examples fun.ProcessSeurat(Var.QuerySeurat)
fun.ProcessSeurat <- function(SeuratObject,nPC=30,Verbose=FALSE){
  SeuratObject %>%
    NormalizeData(verbose = Verbose) %>%
    ScaleData(verbose = Verbose) %>%
    RunPCA(features = rownames(.),npcs=min(nrow(.),100),verbose = Verbose) %>%
    RunUMAP(dims=1:nPC,verbose=Verbose)
}
