brca.se <- base::readRDS(here::here("data_raw","TCGA_SummarizedExperiment_HTseq_BRCA.rds"))

brca.se.data.temp <- brca.se[1:2700, ]
#library(SummarizedExperiment)
gene.annot <-  as.data.frame(SummarizedExperiment::rowData(brca.se.data.temp))
#table(gene.annot$gene_type.)
sample.info <-  as.data.frame(SummarizedExperiment::colData(brca.se.data.temp))

#table(sample.info$year_mda, sample.info$plate_RNAseq)

raw.count <- as.data.frame(SummarizedExperiment::assay(brca.se, 'HTseq_counts'))
fpkm <- as.data.frame(SummarizedExperiment::assay(brca.se, 'HTseq_FPKM'))
fpkm.uq <- as.data.frame(SummarizedExperiment::assay(brca.se, 'HTseq_FPKM.UQ'))

#usethis::use_data(gene.annot, compress = "xz", overwrite = TRUE)
#usethis::use_data(sample.info, compress = "xz", overwrite = TRUE)
#usethis::use_data(raw.count, compress = "xz", overwrite = TRUE)
#usethis::use_data(fpkm, compress = "xz", overwrite = TRUE)
#usethis::use_data(fpkm.uq, compress = "xz", overwrite = TRUE)
usethis::use_data(brca.se.data.temp, compress = "xz", overwrite = TRUE)
