
#' Title
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param ExpreMatrix A matrix of expression levels corresponding to the cell phylogeny.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param Threads Number of cores. (default: 4)
#' @param Redo Redo LUGs searching. (default: FALSE)
#'
#' @returns NULL.
#' @export
#'
#' @examples fun.RunLUEEstimator(Data.TarCA,Var.TPM)
fun.RunLUGEstimator <- function(TarCAObject,ExpreMatrix,DirOut="LEAPOutput/",TopGene=0.1,Threads=4,Redo=FALSE){
  tmp.dirout <- paste0(DirOut,"/LUG/")
  if(!dir.exists(tmp.dirout)){ dir.create(tmp.dirout,recursive = T) }

  tmp.expMatrix <-
  ExpreMatrix %>% t %>%
  .[,which(colMeans(.>0)>TopGene)]

  tmp.tree <- TarCAObject$Tree
  tmp.drop <-
  colnames(tmp.expMatrix) %>%
  mclapply(mc.cores = Threads,function(ggg){
    tmp.fileout <- paste0(tmp.dirout,"/",ggg,".rds")
    if(!file.exists(tmp.fileout)|Redo){
      tmp.ann <-
        data.frame(
          TipLabel=tmp.expMatrix %>% rownames,
          TipAnn=percent_rank(tmp.expMatrix[,ggg])>(1-TopGene)
        )
      tmp.res <- LEU_Estimator(tmp.tree,tmp.ann,ReturnExTree = T)

      if(is.null(tmp.res)){
        tmp.res <- list(Tree=tmp.tree,Ann=tmp.ann,"LEU_Estimator"=data.frame(TipAnn=TRUE,MonoInfo=NA,Total=tmp.ann$TipAnn %>% sum,Np=Inf))
      }

      saveRDS(tmp.res,file=tmp.fileout)
      return()
    }
  })

  message("## Running LUG detector completed.")
}
