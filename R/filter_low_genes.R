# Removing lowly expressed genes

filter.low.genes <- function(data,gene_count,sample_size){
  keep.high <- apply(
    data,
    1,
    function(x) length(x[x>gene_count])>=sample_size
  )
  data.filtered2 <- data[keep.high , ]
  return(data.filtered2)
}
