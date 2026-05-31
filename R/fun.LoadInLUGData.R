#' Title
#'
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param Threads Number of cores. (default: 4)
#'
#' @returns A dataframe of observed LUGs from local files.
#' @export
#'
#' @examples Var.Obs <- fun.LoadInLUGData()
fun.LoadInLUGData <- function(DirOut="LEAPOutput/",Threads=4){
  tmp.ObsLUG <-
    dir(paste0(DirOut,"/LUG/")) %>%
    mclapply(mc.cores = Threads,function(ggg){
      message(paste0(DirOut,"/LUG/",ggg))
      readRDS(paste0(DirOut,"/LUG/",ggg)) %>% .$LEU_Estimator %>% mutate(Symbol=gsub(".rds","",ggg))
    }) %>% bind_rows

  return(tmp.ObsLUG)
}
