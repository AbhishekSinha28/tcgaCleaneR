# RUV-III - PRPS((Pseudo replicate of pseudo sample)) generation

#' @title Generate PRPS for RUV-III
#'
#' @description This function is a part of the data analysis functionality of \code{tgcapkg}. It creates pseudo-replicates
#' of pseudo-samples (PRPS) for unwanted variations like Library size, Batches and Purity in TCGA Pan Cancer Datasets with
#' Cancer biology like Breast Cancer data (BRCA), Lung Cancer (LUAD), Colon Cancer (COAD) and Rectum Cancer (READ). In the
#' function batch refers to the source of Batch Effect variation like Time and Plate which captures variation across
#' biology while factors like Purity captures variation within biology.
#'
#' @param data S4 data object
#' @param librarySize character: Library Size variable in input data object. Default is 'ls'.
#' @param batch character: Batch effect factors. In current package version batch can take values like 'Year', 'Plate' or both
#' @param biology character: Biology of cancer type. TCGA datasets have biology for only four Cancer types i.e.
#' Lung (LUAD), Breast (BRCA), Rectum (READ) & Colon (COAD). So the function supports only these four datasets for RUV-III
#' and PRPS analysis. Default is 'Subtypes'.
#' @param purity character: Purity variable in input data object
#' @param include.ls logical: Do we need to consider library size in creating pseudo samples
#' @param include.purity logical: Do we need to consider purity in creating pseudo samples
#' @param minSamplesPerBatchPS numeric: Minimum number of samples per batch for creating Pseudo Samples
#' @param minSamplesForPuirtyPS numeric: Minimum number of samples for creating Pseudo Samples for purity.
#' @param minSamplesForPurityPerBiology numeric: Number of samples for purity per biology for creating Pseudo Samples
#' @param minSamplesForLibrarySizePerBatch numeric: Number of samples for library size per batch for creating Pseudo Samples
#' @param minSamplesForLibrarySizePS numeric: Minimum number of samples for creating Pseudo Samples for library size
#'
#' @return A S4 list object with the Pseudo replicate for pseudo samples for different batches, library size and purity.
#' @export
#'
#' @examples
#' \dontrun{
#' createPRPS(data=brca.data, batch=c('Year', 'Plates'), purity='Purity_singscore',include.ls=TRUE, include.purity=TRUE,
#' minSamplesPerBatchPS = 3, minSamplesForPuirtyPS = 3, minSamplesForPurityPerBiology = 12,
#' minSamplesForLibrarySizePerBatch = 6,minSamplesForLibrarySizePS = 3)
#' }
createPRPS <- function(data, librarySize = 'ls', batch, biology = 'Subtypes', purity, include.ls, include.purity,
                       minSamplesPerBatchPS, minSamplesForPuirtyPS, minSamplesForPurityPerBiology,
                       minSamplesForLibrarySizePerBatch,minSamplesForLibrarySizePS){

  if(minSamplesForPuirtyPS > minSamplesForPurityPerBiology){
    stop('error: minSamplesForPuirtyPS can not be smaller than minSamplesForPurityPerBiology')
  } else if(minSamplesForPurityPerBiology < 2*minSamplesForPuirtyPS){
    stop('error: minSamplesForPurityPerBiology should be at least two times larger than minSamplesForPuirtyPS')
  } else if(minSamplesForLibrarySizePS > minSamplesForLibrarySizePerBatch) {
    stop('error: minSamplesForLibrarySizePerBatch can not be smaller than minSamplesForLibrarySizePS')
  } else if(minSamplesForLibrarySizePerBatch < 2*minSamplesForLibrarySizePS){
    stop('error: minSamplesForLibrarySizePerBatch should be at least two times larger than minSamplesForLibrarySizePS')
  }

  sample.info <- as.data.frame(SummarizedExperiment::colData(data))
  biology.batch <- c(biology, batch )
  expr.data <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_counts'))
  sample.info[librarySize] <- colSums(expr.data)

  ### Biology
  sample.info$Subtypes <- apply(
    sample.info[, biology, drop = FALSE],
    1,
    paste,
    collapse = "-")

  ### Biology - Batch
  sample.info$biology.batch <- apply(
    sample.info[, biology.batch],
    1, paste,
    collapse = "-")

  ### Create PS per batch
  selected.biology <- unlist(lapply(
    unique(sample.info$Subtypes),
    function(x){
      index <- sample.info$Subtypes == x
      if(sum( table(sample.info$biology.batch[index] ) >= minSamplesPerBatchPS) > 1 ){
        x
      }
    }))
  if(length(selected.biology) > 0){
    message('PRPS are generated for batch effects')
  }else{
    stop('error: there are not enough samples to make pseudo-samples for batch effects, you may want to lower minSamplesPerBatchPS')
  }
  sample.info <- sample.info[sample.info$Subtypes %in% selected.biology , ]
  expr.data <- expr.data[, row.names(sample.info)]

  selected.batches <- names(which(table(sample.info$biology.batch) >= minSamplesPerBatchPS))
  ps.batch <- sapply(
    selected.batches,
    function(x) {
      index <- sample.info$biology.batch == x
      Matrix::rowMeans(expr.data[, index])
    })

  if(include.ls){
    selected.batches.ls <- names(
      which(table(sample.info$biology.batch) >= minSamplesForLibrarySizePerBatch)
    )
    if(length(selected.batches.ls) > 1){
      message('PRPS are generated for library size')
      sample.info <- sample.info[
        with(sample.info,
             order(sample.info[, 'biology.batch'],
                   sample.info[, librarySize])), ]
      expr.data <- expr.data[, row.names(sample.info)]
      ps.ls <- lapply(
        selected.batches.ls,
        function(x){
          index <- sample.info$biology.batch == x
          ls.data <- expr.data[ , index]
          low.ls <- Matrix::rowMeans(ls.data[ , 1:minSamplesForLibrarySizePS])
          high.ls <- rowMeans(ls.data[ , c(ncol(ls.data)-(minSamplesForLibrarySizePS - 1)):ncol(ls.data) ])
          all <- cbind(low.ls, high.ls)
          colnames(all) <- rep(paste(x, 'LS', sep = '-'), 2)
          all
        })
      ps.ls <- do.call(cbind, ps.ls)

    }else{
      stop('error: there are not enough samples to make pseudo-samples for library size effects, you may want to lower minSamplesForLibrarySizePerBatch')
    }
  }else if (!include.ls){
    warning('PRPS is not considered for librray size effects')
    ps.ls = list()
  }

  if(include.purity ){
    selected.biology.purity <- names(
      which(table(sample.info$Subtypes) >= minSamplesForPurityPerBiology)
    )
    if(length(include.purity) > 0){
      message('PRPS are generated for purity effects')
      sample.info <- sample.info[
        with(sample.info,
             order(sample.info[, 'biology.batch'],
                   sample.info[, purity])),]
      expr.data <- expr.data[, row.names(sample.info)]
      ps.purity <- lapply(
        selected.biology.purity,
        function(x) {
          index <- sample.info$Subtypes == x
          purity.data <- expr.data[, index]
          low.pur <- rowMeans(purity.data[, 1:minSamplesForPuirtyPS])
          high.pur <- rowMeans(purity.data[, c(ncol(purity.data) - (minSamplesForPuirtyPS - 1)):ncol(purity.data)])
          all <- cbind(low.pur, high.pur)
          colnames(all) <- rep(paste(x, 'purity', sep = '-'), 2)
          all
        })
      ps.purity <- do.call(cbind, ps.purity)
    }else{
      stop('error: there are not enough samples to make pseudo-samples for purity variation, you may want to lower minSamplesForPurityPerBiology')
    }
  } else if (!include.purity){
    warning('PRPS is not considered for purity effects')
    ps.purity = list()
  }
  return(list(ps.purity = ps.purity, ps.ls = ps.ls, ps.batch = ps.batch))
}
