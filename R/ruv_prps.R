# RUV-III - PRPS((Pseudo replicate of pseudo sample)) generation

#' @title Generate RUV-III PRPS data
#'
#' @description This function is a part of the data analysis functionality of \code{tgcapkg}. It creates pseudo samples for on Library size, Batches and Purity.
#'
#' @param data S4 data object
#' @param batch character: Batch effect factors. In current package version batch can be given values like 'Year', 'Plate' or both.
#' @param biology character: Biology of cancer type. In current package version biology for cancer type Breast, Rectum & Colon is considered.
#' @param purity character: Purity column name
#' @param include.ls logical: Do we need to consider library size in creating pseudo samples.
#' @param include.purity logical: Do we need to consider purity in creating pseudo samples.
#' @param n.ls.ps numeric: Number of samples for library size considered within plate.
#' @param n.sample.ls numeric: Minimum number of samples within library size pseudo sample size.
#' @param n.sample.batch numeric: Number of samples per batch.
#' @param n.sample.purity numeric: Number of samples for purity.
#'
#' @return A S4 list object with the Pseudo replicate for pseudo samples for different batches.
#' @export
#'
#' @examples
#' \dontrun{
#' get.prps(data=brca.data, batch=c('Year', 'Plates'), biology='biology', purity='Purity_singscore', include.ls=TRUE, include.purity=TRUE, n.ls.ps=10, n.sample.batch=3, n.sample.purity=3, n.sample.ls=3)
#' get.prps(data=brca.data, batch=c('Year', 'Plates'), biology='biology', purity=NULL, include.ls=TRUE, include.purity=FALSE, n.ls.ps=10, n.sample.batch=3, n.sample.purity=0, n.sample.ls=3)
#' }

get.prps <- function(data, batch, biology, purity, include.ls, include.purity, n.ls.ps, n.sample.ls,
                     n.sample.batch, n.sample.purity){
  data$log.ls <- log2(colSums(SummarizedExperiment::assay(data, 'HTseq_counts')))
  sample.info <- as.data.frame(SummarizedExperiment::colData(data))
  sample.info$Purity_singscore <- sample.info$purity_HTseq_FPKM
  sample.info$Year <- sample.info$year_mda
  sample.info$Plates <- sample.info$PlateId_mda
  sample.info$TSS <- sample.info$TSS_mda
  sample.info$Tissues <- sample.info$tissue
  sample.info$Center <- sample.info$center_mda
  cols <- c(
    'Year',
    'Plates',
    'TSS',
    'Tissues',
    'Center',
    'log.ls',
    'Purity_singscore'
  )
  sample.info <- sample.info[ , cols]
  sample.info$biology <- sample(letters[1:4], nrow(sample.info), replace = TRUE)
  sample.info$new.batch <- paste0(
    sample.info$Year,
    '_',
    sample.info$Plates
  )
  raw.data <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_counts'))
  gene.annot <- as.data.frame(SummarizedExperiment::rowData(data))
  keep.genes <- gene.annot$gene_biotype_BioMart == 'protein_coding'
  gene.annot <- gene.annot[ keep.genes, ]
  raw.data <- raw.data[ keep.genes, ]
  expr.data = raw.data

  ### combine biology and batch
  sample.info$lib.size <- colSums(expr.data)
  biology.batches <- c(biology, batch )

  ### Biology
  sample.info$biology <- apply(
    sample.info[, biology, drop = FALSE],
    1,
    paste,
    collapse = "-")

  ### Biology - Batch
  sample.info$biology.batches <- apply(
    sample.info[, biology.batches],
    1, paste,
    collapse = "-")

  ### sort samples
  sample.info <- dplyr::arrange(
    sample.info,
    biology.batches,
    lib.size
  )
  expr.data <- expr.data[ , row.names(sample.info)]

  ### Create ps per batch
  selected.batches <-
    names(which(table(sample.info$biology.batches) >= n.sample.batch))
  ps.batch <- sapply(
    selected.batches,
    function(x) {
      index <- sample.info$biology.batches == x
      rowMeans(expr.data[, index])
    })

  if(include.ls){
    selected.batches <- names(
      which(table(sample.info$biology.batches) >= n.ls.ps)
    )
    ps.ls <- lapply(
      selected.batches,
      function(x){
        index <- sample.info$biology.batches == x
        ls.data <- expr.data[ , index]
        low.ls <- rowMeans(ls.data[ , 1:n.sample.ls])
        high.ls <- rowMeans(ls.data[ , c(ncol(ls.data)-(n.sample.ls - 1)):ncol(ls.data) ])
        all <- cbind(low.ls, high.ls)
        colnames(all) <- rep(paste(x, 'LS', sep = '-'), 2)
        all
      })
    ps.ls <- do.call(cbind, ps.ls)
  }else if (! include.ls){
    warning('Library size has not been considered')
    ps.ls = list()
  }

  if(include.purity ){
    selected.biology <- names(
      which(table(sample.info$biology) >= 2*n.sample.purity)
    )
    sample.info <- dplyr::arrange(
      sample.info,
      biology,
      sample.info[ , purity]
    )
    expr.data <- expr.data[ , row.names(sample.info)]
    ps.purity <- lapply(
      selected.biology,
      function(x) {
        index <- sample.info$biology == x
        purity.data <- expr.data[, index]
        low.pur <- rowMeans(purity.data[, 1:n.sample.purity])
        high.pur <- rowMeans(
          purity.data[, c(ncol(purity.data) - (n.sample.purity - 1)):ncol(purity.data)]
        )
        all <- cbind(low.pur, high.pur)
        colnames(all) <- rep(paste(x, 'purity', sep = '-'), 2)
        all

      })
    ps.purity <- do.call(cbind, ps.purity)
  } else if (!include.purity){
    warning('PS do not cover purity variation')
    ps.purity = list()
  }
  return(list(ps.purity = ps.purity, ps.ls = ps.ls, ps.batch = ps.batch))
}


