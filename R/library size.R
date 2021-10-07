# determining library size

library.size <- function(data1,data2,plot_type){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data1))
  library_size <- log2(colSums(data2))
  if (plot_type == "Boxplot"){
    boxplot(library_size ~ sample.info$year_mda)
  } else
    if (plot_type == "Scatterplot"){
      plot(library_size, xlab = 'sample', ylab = 'log2 library size', col = factor(sample.info$year_mda))
    }
}

# df3 <- library.size(data1 = brca.se.data.temp, data2 = df2, plot_type = "Scatterplot")
