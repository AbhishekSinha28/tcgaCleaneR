# Filtering Genes Function

filter.genes <- function(data,gene.type){
  #data.filtered <- vector("logical", length(gene.type))
  data.filtered <- list()
  raw.count <- NULL
  gene.annot.rm <-  as.data.frame(SummarizedExperiment::rowData(data))
  for (i in gene.type){
    keep.genes.rm <- gene.annot.rm$gene_biotype_BioMart == i
    #for (j in length(gene.type)){
    #  data.filtered[j] <- keep.genes.rm
    #}
    data.filtered[[i]] <- data[keep.genes.rm , ]
    #data.filtered[[i]] <- keep.genes.rm
  }
  #data.filtered <- data[keep.genes.rm , ]
  #raw.count.f <- as.data.frame(SummarizedExperiment::assay(data.filtered, 'HTseq_counts'))
  #big_data = dplyr::bind_rows(rbind, data.filtered)
  for (j in 1:length(data.filtered)){
    raw.count.f <- as.data.frame(SummarizedExperiment::assay(data.filtered[[j]], 'HTseq_counts'))
    raw.count <- rbind(raw.count,raw.count.f)
  }
  return(raw.count)
}

#df1 <- rem_genes(data=brca.se.data.temp,gene.type=c("protein_coding","lncRNA"))
