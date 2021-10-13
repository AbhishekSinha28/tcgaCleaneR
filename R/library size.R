# determining library size

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

# df3 <- library.size(data = df2, plot_type = "Scatterplot")
