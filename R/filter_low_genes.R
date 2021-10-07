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

#df2 <- filter.low.genes(data=df1,gene_count = 20,sample_size = 200)
