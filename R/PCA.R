# Generate PCA functioon

get.pca <- function(data, nPcs, is.log){
  .pca <- function(data, nPcs, is.log) {
    if(is.log){
      data <- data
    }else{
      data <- log2(data + 1)
    }
    svd <- BiocSingular::runSVD(
      x = t(data),
      k = nPcs,
      BSPARAM = BiocSingular::bsparam(),
      center = TRUE,
      scale = FALSE
    )
    percent <- svd$d^2/sum(svd$d^2)*100
    percent <-
      sapply(
        seq_along(percent),
        function(i) {round(percent[i], 1)})
    return(list(
      sing.val = svd,
      variation = percent))
  }
  tcga.harmonized <- names(SummarizedExperiment::assays(data))
  pca.cancer.tcga  <- lapply(
    tcga.harmonized,
    function(x){
      .pca(
        data = as.matrix(SummarizedExperiment::assay(data, x)),
        nPcs = nPcs,
        is.log = FALSE)
    })
  names(pca.cancer.tcga) <- tcga.harmonized
  return(pca.cancer.tcga)
}

#df6 <- get.pca(data = df5, nPcs = 10, is.log = FALSE)
#df7 <- get.pca(data = df5, nPcs = 8, is.log = FALSE)
