# determining library size

#' @title Library Size
#'
#' @description This function is a part of the data analysis functionality of tgcapkg.    It allows user to analyze the Library Size bias, a technical bias in the TCGA RNA-seq.    The user can input the TGCA dataset and the type of plot to analyse library size.
#'
#' @param data Input TGCA Dataset.
#' @param plot_type character: Plot type
#'
#' @return Scatter Plot, Box plot
#' @export
#'
#' @examples
#' library.size(data = brca.data, plot_type = "Scatterplot")
#' \dontrun{
#' library.size(data = brca.data, plot_type = "Boxplot")
#' }
library.size <- function(data,plot_type){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library_size <- log2(colSums(raw.count))
  if (plot_type == "Boxplot"){
    boxplot(library_size ~ sample.info$year_mda)
  } else
    if (plot_type == "Scatterplot"){
      plot(library_size, xlab = 'sample', ylab = 'log2 library size', col = factor(sample.info$year_mda))
    }
}
