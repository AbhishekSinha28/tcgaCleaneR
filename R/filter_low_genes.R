# Removing lowly expressed genes

#' @title Low Count Genes Filter
#'
#' @description This function is a part of the data wrangling functionality of tgcapkg.    It allows user to input the TGCA dataset and the threshold for the minimum gene count and sample count.
#'
#' @param data Input TGCA Dataset.
#' @param gene_count numeric: gene count threshold
#' @param sample_size numeric: sample size threshold
#'
#' @return S4 data object
#' @export
#'
#' @examples
#' filter.low.genes(data=brca.data,gene_count = 20,sample_size = 200)
filter.low.genes <- function(data,gene_count,sample_size){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  keep.high <- apply(
    raw.count,
    1,
    function(x) length(x[x>gene_count])>=sample_size
  )
  data.filtered <- data[keep.high , ]
  return(data.filtered)
}
