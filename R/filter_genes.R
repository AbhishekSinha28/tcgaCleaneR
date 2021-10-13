# Filtering Genes Function

filter.genes <- function(data,gene.type){
  gene.annot.rm <-  as.data.frame(SummarizedExperiment::rowData(data))
  keep.genes.rm <- gene.annot.rm$gene_biotype_BioMart %in% gene.type
  data.filtered <- data[keep.genes.rm , ]
  return(data.filtered)
}

#df1 <- filter.genes(data=brca.se.data.temp,gene.type=c("protein_coding","lncRNA"))
