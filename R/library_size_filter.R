# Removing samples based on library size

library.size.fil <- function(data,ls_cutoff){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library_size <- log2(colSums(raw.count))
  keep.samples <- library_size > ls_cutoff
  brca.se.filtered <- data[ , keep.samples]
  return(brca.se.filtered)
}

# df4 <- library.size.fil(data = df2, ls_cutoff = 22)
