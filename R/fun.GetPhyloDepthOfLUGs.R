#' Title
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param LUG A dataframe of LUGs obtained from fun.ExtractLUG.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param Threads Number of cores. (default: 4)
#'
#' @returns A dataframe of normalized depth of internal nodes with LUG activation events.
#' @export
#'
#' @examples fun.GetPhyloDepthOfLUGs(Data.TarCA,Var.LUGs)
fun.GetPhyloDepthOfLUGs <- function(TarCAObject,LUG,DirOut="LEAPOutput/",Threads=4){
  Var.PhyloDepthOfLUGs <-
  LUG %>% filter(Sig) %>% .$Symbol %>%
  mclapply(mc.cores = Threads,function(ggg){
    tmp2 <- TarCAObject$Name2Meta %>% mutate(Nor_x=x/median(x[isTip]))

    tmp1 <- readRDS(paste0(DirOut,"/LUG/",ggg,".rds"))
    tmp3 <-
    tmp1$FilterBiasParent %>%
    filter(Filter) %>%
    rownames %>%
    lapply(function(nnn){
      tmp1$AllDescendants[[nnn]] %>%
        filter(TipAnn) %>% select(TipLabel) %>%
        left_join(TarCAObject$Ann) %>%
        mutate(Node=nnn,Nor_x=tmp2[Node,"Nor_x"])
    }) %>% bind_rows %>%
    mutate(Symbol=ggg)
    return(tmp3)
  }) %>% bind_rows

  return(Var.PhyloDepthOfLUGs)
}
