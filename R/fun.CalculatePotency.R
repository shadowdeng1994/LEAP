#' @title Calculating potency for each internal node..
#' @description
#' Mapping inferred state to the nodes and their unsolved descendants.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Nodes A vector of internal nodes.
#' @param Threads Number of cores. (default: 4)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#'
#' @returns A dataframe with potency information for each internal node.
#' @export
#'
#' @examples Var.Potency <- CalculatePotency(Data.TarCA,Var.PredictedState$NodeLabel)
fun.CalculatePotency <- function(TarCAObject,Nodes,Threads=4,DirOut="LEAPOutput/",SaveObject=TRUE){
  Var.Potency <-
  Nodes %>%
  mclapply(mc.cores = Threads,function(nnn){
    tmp1 <- TarCAObject$AllDescendants[[nnn]]$TipAnn %>% sort %>% as.character %>% table

    tmp2 <- tmp1/sum(tmp1)
    tmp.SEntropy <- -sum(tmp2*log2(tmp2))

    data.frame(
      NodeName=nnn,
      CladeSize=sum(tmp1),
      On=length(tmp1),
      ShannonEntropy=tmp.SEntropy,
      DetailInfo=paste(paste0(names(tmp1),"[",tmp1,"]"),collapse = ";")
    )
  }) %>% bind_rows

  if(SaveObject){
    fun.message("==> Saving Potency to ",DirOut,".")
    saveRDS(Var.Potency,paste0(DirOut,"/Var.Potency.rds"))
  }

  return(Var.Potency)
}
