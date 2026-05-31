#' @title Run label transfer from reference data to query data.
#' @description
#' Run label transfer from reference data to query data, both of which were binarized by LUG.
#'
#' @param QueryMatrix A binarized matrix of query data obtained from GetInternalMatrix.
#' @param RefMatrix A binarized matrix of reference data obtained from GetReferenceMatrix.
#' @param TransferLabel Transfer Labels from reference to query.
#' @param QuerySeurat Existed QuerySeurat.
#' @param RefSeurat Existed RefSeurat.
#' @param Anchors Existed Anchor.
#' @param nPC Number of PC. (default: 30)
#' @param Verbose Verbose or not. (default: FALSE)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#'
#' @returns A dataframe of transferred result.
#' @export
#'
#' @examples Var.Predictions <- RunLabelTransfer(Var.InternalMatrix,Var.ReferenceMatrix,Data.Ref_metadata$lineage)
RunLabelTransfer <- function(QueryMatrix=NULL,RefMatrix=NULL,TransferLabel,QuerySeurat=NULL,RefSeurat=NULL,Anchors=NULL,nPC=30,Verbose=FALSE,DirOut="LEAPOutput/",SaveObject=TRUE){
  if(is.null(Anchors)){
    if(is.null(QuerySeurat)|is.null(RefSeurat)){
      fun.message("--> Creating Seurat object.")
      Var.QuerySeurat <- CreateSeuratObject(counts=QueryMatrix,project = "Query")
      Var.RefSeurat <- CreateSeuratObject(counts=RefMatrix,project = "Ref")
      # ,meta.data = RefMetadata

      fun.message("==> Processing Seurat object.")
      Var.QuerySeurat <- Var.QuerySeurat %>% fun.ProcessSeurat(nPC,Verbose)
      Var.RefSeurat <- Var.RefSeurat %>% fun.ProcessSeurat(nPC,Verbose)

      if(SaveObject){
        fun.message("==> Saving Seurat to ",DirOut,".")
        saveRDS(Var.QuerySeurat,paste0(DirOut,"/Var.QuerySeurat.rds"))
        saveRDS(Var.RefSeurat,paste0(DirOut,"/Var.RefSeurat.rds"))
      }
    }else{
      fun.message("==> Loading Seurat object.")
      Var.QuerySeurat <- QuerySeurat
      Var.RefSeurat <- RefSeurat
    }

    fun.message("==> Finding transfer anchors.")
    Var.Anchors <- FindTransferAnchors(reference = Var.RefSeurat, query = Var.QuerySeurat, dims = 1:nPC,reference.reduction = "pca",features = rownames(Var.QuerySeurat),verbose=Verbose)

    if(SaveObject){
      fun.message("==> Saving Anchor to ",DirOut,".")
      saveRDS(Var.Anchors,paste0(DirOut,"/Var.Anchors.rds"))
    }
  }else{
    fun.message("==> Loading transfer anchors.")
    Var.Anchors <- Anchors
  }

  fun.message("==> Running label transfer.")
  Var.Predictions <- TransferData(anchorset = Var.Anchors, refdata = TransferLabel, dims = 1:nPC,verbose=Verbose)

  if(SaveObject){
    fun.message("==> Saving Prediction to ",DirOut,".")
    saveRDS(Var.Predictions,paste0(DirOut,"/Var.Predictions.rds"))
  }

  return(Var.Predictions)
}
