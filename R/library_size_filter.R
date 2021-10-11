# Removing samples based on library size

library.size.fil <- function(data,ls_cutoff){
  library_size <- log2(colSums(data))
  keep.samples <- library_size > ls_cutoff
  brca.se.filtered3 <- data[ , keep.samples]
  return(brca.se.filtered3)
}

# df4 <- library.size.fil(data = df2, ls_cutoff = 22)
