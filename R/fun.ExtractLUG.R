#' @title Extract the LUGs from cell phylogeny and corresponding transcritptome.
#' @description
#' LUGs inference based on permutation test.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param ExpreMatrix A matrix of expression levels corresponding to the cell phylogeny.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param NumShuf Number of shuffling. (default: 1000)
#' @param p_threshold Threshold of empirical P-value. (default: 0.01)
#' @param SeedFrom Set the beginning of seed. (default: NULL)
#' @param Threads Number of cores. (default: 4)
#' @param Redo Redo LUGs searching. (default: FALSE)
#'
#' @returns A dataframe of inferred LUGs.
#' @export
#'
#' @examples Data.TarCA <- CreateTarCAObject(Var.Tree,Var.Ann)
#' @examples Var.LUGs <- fun.ExtractLUG(Data.TarCA,Var.TPM)
fun.ExtractLUG <- function(TarCAObject,ExpreMatrix,DirOut="LEAPOutput/",TopGene=0.1,NumShuf=1000,p_threshold=0.01,SeedFrom=NULL,Threads=4,Redo=FALSE){
  if(!dir.exists(DirOut)){
    fun.message("--> ",DirOut,"is absent. Creating.")
    dir.create(DirOut)
  }

  fun.message("--> Running LUG estimator.")
  fun.RunLUGEstimator(TarCAObject,ExpreMatrix,DirOut,TopGene,Threads,Redo)
  Var.Obs <- fun.LoadInLUGData(DirOut,Threads)

  fun.message("--> Shuffling (",NumShuf,"times).")
  Var.Shuf <- fun.RunShufLUGEstimator(TarCAObject,Var.Obs,DirOut,NumShuf,SeedFrom,Threads,Redo)

  fun.message("--> Extracting high confident LUGs.")
  Var.LUGs <- fun.DetLUG(Var.Obs,Var.Shuf,p_threshold,Threads)

  return(Var.LUGs)
}
