# Purity filter function

purity.filter <- function(data,purity_cutoff){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  keep.samples <- sample.info$purity_HTseq_counts > purity_cutoff
  brca.se.filtered <- data[keep.samples,]
  return(brca.se.filtered)
}

#df3 <- purity.filter(data=df2,purity_cutoff= 0.496)
