#' Title
#'
#' @param PhyloDepth A dataframe of PhyloDepth obtained from fun.GetPhyloDepthOfLUGs.
#' @param MinAct Minimum number of activation events for reconstructed internal nodes. (default: 3)
#' @param Threads Number of cores. (default: 4)
#'
#' @returns A matrix of inferred transcriptome of internal nodes.
#' @export
#'
#' @examples fun.GetInternalState(Var.PhyloDepth)
fun.GetInternalState <- function(PhyloDepth,MinAct=3,Threads=4){
  tmp1 <- PhyloDepth %>% split(x=.,f=.$TipLabel)

  tmp4 <-
    names(tmp1) %>%
    mclapply(mc.cores = Threads,function(nnn){
      tmp1[[nnn]] %>%
        arrange(TipLabel,Nor_x) %>%
        select(TipLabel,Node,Symbol,Node_x=Nor_x,Symbol_x=Nor_x) %>%
        spread(Symbol,Symbol_x) %>% gather(Symbol,Symbol_x,-TipLabel,-Node,-Node_x) %>%
        group_by(TipLabel,Symbol) %>% mutate(Symbol_x=replace_na(Symbol_x,Symbol_x[!is.na(Symbol_x)][1])) %>%
        group_by %>% mutate(isAct=Symbol_x<=Node_x) %>%
        select(Node,Symbol,isAct)
    }) %>% bind_rows %>%
    unique

  Var.InternalState <-
    tmp4 %>%
    group_by(Node) %>% filter(sum(isAct)>=MinAct) %>%
    group_by %>%
    mutate(isAct=as.numeric(isAct)) %>%
    spread(Symbol,isAct,0) %>%
    column_to_rownames("Node") %>%
    t

  return(Var.InternalState)
}
