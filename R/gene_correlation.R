#  Gene Correlation analysis function

gene.corr <- function(data, is.log, type, cor.method, n.cores){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
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
  }
  sample.info$ls <- log2(colSums(raw.count))
  raw.count <- as.matrix(raw.count)
  if(type == "librarysize"){
    genes.df <- .corr.gene.variable(
      expr.data = raw.count,
      is.log = is.log,
      variable = sample.info$ls,
      method = cor.method,
      n.cores = 5
    )
    hist(genes.df$rho)
  } else
    if(type == "purity_HTseq_counts"){
      genes.df <- .corr.gene.variable(
        expr.data = raw.count,
        is.log = is.log,
        variable = sample.info$purity_HTseq_counts,
        method = cor.method,
        n.cores = 5
      )
      hist(genes.df$rho)
    } else
      if(type == "purity_HTseq_FPKM"){
        genes.df <- .corr.gene.variable(
          expr.data = raw.count,
          is.log = is.log,
          variable = sample.info$purity_HTseq_FPKM,
          method = cor.method,
          n.cores = 5
        )
        hist(genes.df$rho)
      } else
        if(type == "purity_HTseq_FPKM.UQ"){
          genes.df <- .corr.gene.variable(
            expr.data = raw.count,
            is.log = is.log,
            variable = sample.info$purity_HTseq_FPKM.UQ,
            method = cor.method,
            n.cores = 5
          )
          hist(genes.df$rho)
        }
  return(genes.df)
}

#df9 <- gene.corr(data = df5, is.log = FALSE, type = "librarysize", cor.method = 'spearman', n.cores = 5)
#gene.corr(data = df5, is.log = FALSE, type = "purity_HTseq_FPKM", cor.method = 'spearman', n.cores = 5)
#gene.corr(data = df5, is.log = FALSE, type = "purity_HTseq_FPKM", cor.method = 'pearson', n.cores = 5)
