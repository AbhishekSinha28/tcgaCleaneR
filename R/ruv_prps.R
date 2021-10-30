.CreatePseudoSamplesForLsPurityBatch <- function(
  expr.data,
  sample.info,
  batch,
  biology,
  purity,
  include.ls,
  include.purity,
  n.ls.ps,
  n.sample.batch,
  n.sample.purity,
  n.sample.ls){
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
    collapse = "-"
  )
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
    print('Library size has not been considered')
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
    print('PS do not cover purity variation')
    ps.purity = list()
  }
  return(list(ps.purity = ps.purity, ps.ls = ps.ls, ps.batch = ps.batch))
}

raw.data <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_counts'))
gene.annot <- as.data.frame(SummarizedExperiment::rowData(data))
keep.genes <- gene.annot$gene_biotype_BioMart == 'protein_coding'
gene.annot <- gene.annot[ keep.genes, ]
raw.data <- raw.data[ keep.genes, ]

prps <-
  .CreatePseudoSamplesForLsPurityBatch(
    expr.data = raw.data,
    sample.info = sample.info,
    batch = c('Year', 'Plates'),
    biology = 'biology',
    purity = NULL,
    include.ls = T,
    include.purity = F,
    n.ls.ps = 10,
    n.sample.batch = 3,
    n.sample.purity = 3,
    n.sample.ls = 3
  )
prps$ps.batch
prps.batch <- prps$ps.batch
colnames(prps.batch) <- unlist(lapply(
  colnames(prps.batch),
  function(x) strsplit(x, '-')[[1]][1]
))
