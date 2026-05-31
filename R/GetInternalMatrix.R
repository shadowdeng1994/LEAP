#' @title Reconstruct the LUG-transcriptome of internal nodes.
#' @description
#' Reconstruct the transcriptome of internal nodes based on LUG inference.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject. (default: NULL)
#' @param ExpreMatrix A matrix of expression levels corresponding to the cell phylogeny. (default: NULL)
#' @param LUGs A dataframe of LUGs obtained from fun.ExtractLUG. (optional,default: NULL)
#' @param PhyloDepth A dataframe of PhyloDepth obtained from fun.GetPhyloDepthOfLUGs. (optional,default: NULL)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param NumShuf Number of shuffling. (default: 1000)
#' @param p_threshold Threshold of empirical P-value. (default: 0.01)
#' @param MinAct Minimum number of activation events for reconstructed internal nodes. (default: 3)
#' @param SeedFrom Set the beginning of seed. (default: NULL)
#' @param Threads Number of cores. (default: 4)
#' @param SaveObject Saving all objects. (default: TRUE)
#' @param Redo Redo LUGs searching. (default: FALSE)
#'
#' @returns A matrix of inferred transcriptome of internal nodes.
#' @export
#'
#' @examples Data.TarCA <- CreateTarCAObject(Var.Tree,Var.Ann)
#' @examples Var.InternalState <- GetInternalMatrix(Data.TarCA,Var.TPM)
GetInternalMatrix <- function(TarCAObject=NULL,ExpreMatrix=NULL,LUGs=NULL,PhyloDepth=NULL,DirOut="LEAPOutput/",TopGene=0.1,NumShuf=1000,p_threshold=0.01,MinAct=3,SeedFrom=NULL,Threads=4,SaveObject=TRUE,Redo=FALSE){
  fun.message("==> Extracting LUG.")
  Var.LUGs <- fun.ExtractLUG(TarCAObject,ExpreMatrix,DirOut,TopGene,NumShuf,p_threshold,SeedFrom,Threads,Redo)

  fun.message("==> Converting to PhyloDepth.")
  Var.PhyloDepthOfLUGs <- fun.GetPhyloDepthOfLUGs(TarCAObject,Var.LUGs,DirOut="LEAPOutput/",Threads)

  fun.message("==> Inferring internal state.")
  Var.InternalState <- fun.GetInternalState(Var.PhyloDepthOfLUGs,MinAct,Threads)

  if(SaveObject){
    fun.message("==> Saving to ",DirOut,".")
    saveRDS(Var.LUGs,paste0(DirOut,"/Var.LUGs.rds"))
    saveRDS(Var.PhyloDepthOfLUGs,paste0(DirOut,"/Var.PhyloDepthOfLUGs.rds"))
    saveRDS(Var.InternalState,paste0(DirOut,"/Var.InternalState.rds"))
  }

  return(Var.InternalState)
}


