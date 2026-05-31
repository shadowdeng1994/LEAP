#' Title
#'
#' @param TarCAObjeLUG A TarCA object obtained from CreateTarCAObject.
#' @param ObsLUG A dataframe of LUGs obtained from fun.LoadInLUGData.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param NumShuf Number of shuffling. (default: 1000)
#' @param SeedFrom Set the beginning of seed. (default: NULL)
#' @param Threads Number of cores. (default: 4)
#' @param Redo Redo LUGs searching. (default: FALSE)
#'
#' @returns A dataframe of shuffling LUGs.
#' @export
#'
#' @examples Var.Shuf <- fun.RunShufLUGEstimator(Data.TarCA,Var.Obs)
fun.RunShufLUGEstimator <- function(TarCAObject,ObsLUG,DirOut="LEAPOutput/",NumShuf=1000,SeedFrom=NULL,Threads=4,Redo=FALSE){
  tmp.dirout <- paste0(DirOut,"/Shuffle/")
  if(!dir.exists(tmp.dirout)){ dir.create(tmp.dirout,recursive = T) }

  tmp.tree <- TarCAObject$Tree

  tmp.set <- ObsLUG %>% .$Total %>% unique %>% sort
  tmp.shuf <-
  tmp.set %>%
  lapply(function(sss){
    fun.message("--> Total: ",sss,"(",match(sss,tmp.set),"/",length(tmp.set),")")

    tmp.fileout <- paste0(tmp.dirout,"/Total_",sss,".rds")

    if(!file.exists(tmp.fileout)|Redo){
      message("\tRunning...")

      tmp.out <-
        1:NumShuf %>%
        mclapply(mc.cores = Threads,function(rrr){
          if(!is.null(SeedFrom)){ set.seed(SeedFrom+rrr) }
          tmp.ann <- data.frame(TipLabel=tmp.tree$tip.label) %>% mutate(TipAnn=sample(1:n()<=sss))
          tmp.res <- LEU_Estimator(tmp.tree,tmp.ann)
          if(is.null(tmp.res)){
            tmp.res <- data.frame(TipAnn=TRUE,MonoInfo=NA,Total=tmp.ann$TipAnn %>% sum,Np=Inf)
          }
          return(tmp.res)
        }) %>% bind_rows

      saveRDS(tmp.out,file=tmp.fileout)
    }else{
      message(paste0("\tLoading: ",tmp.fileout))
      tmp.out <- readRDS(tmp.fileout)
    }
    return(tmp.out)
  }) %>% bind_rows

  return(tmp.shuf)
}
