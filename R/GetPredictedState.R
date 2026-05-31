#' @title Get Inferred states from label transfer data.
#' @description
#' Mapping inferred state to the nodes and their unsolved descendants.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Predictions A dataframe of transferred result.
#' @param RemoveRoot Remove root's name.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#'
#' @returns A datafraome of impuated annotation on the cell phylogeny.
#' @export
#'
#' @examples Var.PredictedState <- GetPredictedState(Data.TarCA,Var.Predictions)
GetPredictedState <- function(TarCAObject,Predictions,RemoveRoot=NULL,DirOut="LEAPOutput/",SaveObject=TRUE){
  Var.PredictedState <-
  fun.AssignLaterState(
    TarCAObject,
    Predictions %>% select(TransLabel=predicted.id)
  ) %>%
  mutate(TransLabel=Predictions[AncLabel,"predicted.id"])

  if(!is.null(RemoveRoot)){ Var.PredictedState <- Var.PredictedState %>% filter(TransLabel!="Root") }

  if(SaveObject){
    fun.message("==> Saving PredictedState to ",DirOut,".")
    saveRDS(Var.PredictedState,paste0(DirOut,"/Var.PredictedState.rds"))
  }

  return(Var.PredictedState)
}
