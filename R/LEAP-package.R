#' LEAP: LUG-Encoded Ancestral Projection
#'
#' A phylogeny-based computational framework for inferring the transcriptional
#' states of internal nodes from cell phylogenies via lineage-specific
#' upregulated genes (LUGs).
#'
#' @importFrom TarCA.beta Np_Estimator LEU_Estimator
#' @importFrom dplyr filter mutate select arrange bind_rows left_join group_by
#' @importFrom tidyr spread gather replace_na
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr %>%
#' @importFrom Seurat CreateSeuratObject FindTransferAnchors TransferData NormalizeData ScaleData RunPCA RunUMAP
#' @importFrom parallel mclapply
#' @importFrom stats median
"_PACKAGE"
