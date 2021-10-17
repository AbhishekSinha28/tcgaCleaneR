# Vector correlation analysis function

pca.corr <- function(pca.data, data, type){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  time.years <- fastDummies::dummy_cols(sample.info$year_mda)
  time.years <- time.years[, c(2:ncol(time.years))]
  data.set.names <- names(SummarizedExperiment::assays(data))
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  nPCs <- 10
  cca.time.year <-
    lapply(
      data.set.names,
      function(x){
        sapply(
          1:nPCs,
          function(y) {
            cca.pam50 <- stats::cancor(
              x = pca.data[[x]]$sing.val$u[, 1:y, drop = FALSE],
              y = time.years)
            1 - prod(1 - cca.pam50$cor ^ 2)
          })

      })
  .corr.gene.variable <- function(expr.data, is.log, variable, method, n.cores){
    if(is.log){
      expr.data <- expr.data
    }else{
      expr.data <- log2(expr.data + 1)
    }
    rho <- parallel::mclapply(
      1:nrow(expr.data),
      function(x){
        round(cor.test(x = expr.data[x, ], y = variable, method = method)[[4]], 6)},
      mc.cores = n.cores)

    pval <- parallel::mclapply(
      1:nrow(expr.data),
      function(x){
        cor.test(x = expr.data[x, ], y = variable, method = method)[[3]]},
      mc.cores = n.cores)

    results <- data.frame(
      genes = row.names(expr.data),
      rho = unlist(rho),
      pvalue = unlist(pval),
      adj.pvalue = p.adjust(unlist(pval), 'BH')
    )
    return(results)
  }
  raw.count <- as.matrix(raw.count)
  if(type == "librarysize"){
    sample.info$ls <- log2(colSums(raw.count))
    genes.ls <- .corr.gene.variable(
      expr.data = raw.count,
      is.log = FALSE,
      variable = sample.info$ls,
      method = 'spearman',
      n.cores = 5
    )
    hist(genes.ls$rho)
  } else
    if(type == "purity"){
      genes.purity <- .corr.gene.variable(
        expr.data = raw.count,
        is.log = FALSE,
        variable = sample.info$purity_HTseq_FPKM,
        method = 'spearman',
        n.cores = 5
      )
      hist(genes.purity$rho)
    }

}

#pca.corr(pca.data = df5, data = df4, type = "librarysize")
