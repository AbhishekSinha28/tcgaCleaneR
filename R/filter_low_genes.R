# Removing lowly expressed genes

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

#df2 <- filter.low.genes(data=df1,gene_count = 20,sample_size = 200)
