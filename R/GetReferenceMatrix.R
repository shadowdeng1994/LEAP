#' @title Binarize reference dataset.
#' @description
#' Convert reference dataset into binarized matrix.
#'
#' @param ExpreMatrix A matrix of expression levels corresponding to the reference atlas.
#' @param Metadata A dataframe containing metadata for each reference cell.
#' @param GroupBy Column name used to group reference data. (default: NULL)
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param Threads Number of cores. (default: 4)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#'
#' @returns A binarized matrix of reference data.
#' @export
#'
#' @examples GetReferenceMatrix(Data.Ref_Count,Data.Ref_metadata,GroupBy="embryo.time.bin")
GetReferenceMatrix <- function(ExpreMatrix,Metadata,GroupBy=NULL,TopGene=0.1,Threads=4,DirOut="LEAPOutput/",SaveObject=TRUE){
  if(!is.null(GroupBy)){
    tmp.set <- split(x=rownames(Metadata),f=Metadata[,GroupBy])
  }else{
    tmp.set <- list("All"=rownames(Metadata))
  }

  Var.RefMatrix <-
  names(tmp.set) %>%
  lapply(function(ttt){
    fun.message("--> Group: ",ttt)

    tmp2 <- tmp.set[[ttt]] %>% ExpreMatrix[,.]

    tmp3 <-
    1:nrow(tmp2) %>%
    mclapply(mc.cores = Threads,function(iii){
      as.numeric(percent_rank(as.numeric(tmp2[iii,]))>(1-TopGene))
    }) %>%
    do.call(rbind,.) %>%
    as.data.frame

    dimnames(tmp3) <- dimnames(tmp2)

    return(tmp3)
  }) %>% bind_cols

  if(SaveObject){
    fun.message("==> Saving to ",DirOut,".")
    saveRDS(Var.RefMatrix,paste0(DirOut,"/Var.RefMatrix.rds"))
  }

  return(Var.RefMatrix)
}
