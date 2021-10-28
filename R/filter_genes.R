# Filtering Genes Function

#' @title Gene Filter
#'
#' @description This function is a part of the data wrangling functionality of tgcapkg.    It allows user to input the TGCA dataset and the required genes to filter the input data based on genes.
#'
#' @usage gene.filter(data,gene.type)
#'
#' @param data Input TGCA Dataset.
#' @param gene.type A character vector of items.
#'
#' @return S4 data object
#' @export
#'
#' @examples
#' \dontrun{
#' gene.filter(data=brca.data,gene.type=c("protein_coding","lncRNA"))
#' gene.filter(data=brca.data,gene.type=c("protein_coding"))
#' }
gene.filter <- function(data,gene.type){
  gene.annot.rm <-  as.data.frame(SummarizedExperiment::rowData(data))
  keep.genes.rm <- gene.annot.rm$gene_biotype_BioMart %in% gene.type
  data.filtered <- data[keep.genes.rm , ]
  return(data.filtered)
}

