#' @title Create a TarCA object
#' @description
#' Creating a TarCA object with the tree object and the annotation dataset of tips
#'
#' @param Tree Tree object from read.tree.
#' @param Ann Dataframe with columns TipLabel and TipAnn.
#' @param FileOut Path to saving the TarCA object. (default: NULL)
#'
#' @returns A TarCA object.
#' @export
#'
#' @examples Data.TarCA <- CreateTarCAObject(Var.Tree,Var.Ann)
CreateTarCAObject <- function(Tree,Ann,FileOut=NULL){
  Data.TarCA <- Np_Estimator(Tree,Ann,ReturnExTree = T)
  if(!is.null(FileOut)){
    save(Data.TarCA,file=FileOut)
  }
  return(Data.TarCA)
}
