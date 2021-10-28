# Purity filter function

#' @title Filter Samples based on Tumor Purity
#'
#' @description This function is a part of the data wrangling functionality of tgcapkg.    It allows user to handle sample purity in TCGA Cancer data by filtering out the samples that are above a specific purity threshold.
#'
#' @param data Input TGCA Dataset.
#' @param purity_cutoff numeric: Sample purity cutoff
#'
#' @return S4 data object
#' @export
#'
#' @examples
#' \dontrun{
#' purity.filter(data= brca.data,purity_cutoff= 0.496)
#' }
purity.filter <- function(data,purity_cutoff){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  keep.samples <- sample.info$purity_HTseq_counts > purity_cutoff
  brca.se.filtered <- data[ ,keep.samples]
  return(brca.se.filtered)
}

