#' Title
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Predictions A dataframe of transferred result.
#'
#' @returns A dataframe with later state.
#' @export
#'
#' @examples fun.AssignLaterStage(Data.TarCA,Var.Predictions)
fun.AssignLaterState <- function(TarCAObject,Predictions){
  fun.looping <- function(TarCAObject,Predictions,NNNode,PPParentState,AAAncLUGNode){
    tmp.state <- Predictions[NNNode,"TransLabel"]
    tmp.group <- "LUG"
    tmp.anclabel <- AAAncLUGNode
    if(is.na(tmp.state)){ tmp.state <- PPParentState;tmp.group <- "Legacy" }
    if(tmp.group=="LUG"){ tmp.anclabel <- NNNode }
    tmp.out <- data.frame(NodeLabel=NNNode,TransLabel=tmp.state,Type=tmp.group,AncLabel=tmp.anclabel)

    tmp.daughter <- TarCAObject$Name2Daughter[[NNNode]] %>% TarCAObject$Name2Meta[.,] %>% filter(!isTip) %>% rownames
    if(length(tmp.daughter)>0){
      tmp.out <-
        tmp.daughter %>%
        lapply(function(nnn){
          fun.looping(TarCAObject,Predictions,nnn,tmp.state,tmp.anclabel)
        }) %>% bind_rows %>%
        rbind(tmp.out,.)
    }

    return(tmp.out)
  }

  tmp.res <-
    fun.looping(
      TarCAObject,
      Predictions,
      TarCAObject$Name2Meta %>% filter(isRoot) %>% rownames,
      "Root",
      TarCAObject$Name2Meta %>% filter(isRoot) %>% rownames
    )

  return(tmp.res)
}
