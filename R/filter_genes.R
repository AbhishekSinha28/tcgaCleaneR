# Filtering Genes Function

filter.genes <- function(data,gene.type){
  #data.filtered <- vector("logical", nrow(data))
  gene.annot.rm <-  as.data.frame(SummarizedExperiment::rowData(data))
  for (i in gene.type){
    keep.genes.rm <- gene.annot.rm$gene_biotype_BioMart == i
    #data.filtered <- data[keep.genes.rm , ]
  }
  data.filtered <- data[keep.genes.rm , ]
  raw.count.f <- as.data.frame(SummarizedExperiment::assay(data.filtered, 'HTseq_counts'))
  return(raw.count.f)
}
