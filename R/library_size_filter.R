# Removing samples based on library size

#' @title Filter Samples based on Library Size
#'
#' @description This function is a part of the data wrangling functionality of tgcapkg.    It allows user to handle the bias in TCGA Cancer data due to library size by filtering out the samples with sample size greater than the threshold. Using \code{library.size}, user can determine the threshold.
#'
#' @param data Input TGCA Dataset.
#' @param ls_cutoff numeric: library size threshold
#'
#' @return S4 data object
#' @export
#'
#' @examples
#'
#' library.size.fil(data = brca.data, ls_cutoff = 17.5)
#'
library.size.fil <- function(data,ls_cutoff){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library_size <- log2(colSums(raw.count))
  keep.samples <- library_size > ls_cutoff
  brca.se.filtered <- data[ , keep.samples]
  return(brca.se.filtered)
}


