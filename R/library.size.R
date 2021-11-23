# determining library size

#' @title Library Size
#'
#' @description This function is a part of the data analysis functionality of tgcapkg.
#' It allows user to analyze the Library Size bias, a technical bias in the TCGA RNA-seq.
#' The user can input the \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) and the type of plot
#' to analyse variations in library size across years and samples.
#'
#' @param data SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.
#' @param plot_type character: Plot type
#'
#' @return Scatter Plot, Box plot
#' @export
#'
#' @examples
#' plotLibSize(data = brca.data, plot_type = "Scatterplot")
#' \dontrun{
#'
#' plotLibSize(data = brca.data, plot_type = "Boxplot")
#' }
plotLibSize <- function(data,plot_type){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library_size <- log2(colSums(raw.count))
  if (plot_type == "Boxplot"){
    boxplot(library_size ~ sample.info$Year, xlab = 'Sample Years', ylab = 'log2 library size')
  } else
    if (plot_type == "Scatterplot"){
      plot(library_size, xlab = 'sample', ylab = 'log2 library size', col = factor(sample.info$Year))
      legend("bottomright", legend=levels(factor(sample.info$Year)), fill = 1:5, #as.numeric(levels(factor(sample.info$Year))),
             title="Year",
             inset=c(0,0.92), xpd=TRUE, horiz=TRUE, bty="n"
      )
    }
}
