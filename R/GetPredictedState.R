#' @title Get Inferred states from label transfer data.
#' @description
#' Mapping inferred state to the nodes and their unsolved descendants.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Predictions A dataframe of transferred result.
#' @param RemoveRoot Remove root's name.
#'
#' @returns A datafraome of impuated annotation on the cell phylogeny.
#' @export
#'
#' @examples GetPredictedState(Data.TarCA,Var.Predictions)
GetPredictedState <- function(TarCAObject,Predictions,RemoveRoot=NULL){
  Var.PredictedState <-
  fun.AssignLaterState(
    TarCAObject,
    Predictions %>% select(TransLabel=predicted.id)
  ) %>%
  mutate(TransLabel=Predictions[AncLabel,"predicted.id"])

  if(!is.null(RemoveRoot)){ Var.PredictedState <- Var.PredictedState %>% filter(TransLabel!="Root") }

  return(Var.PredictedState)
}
