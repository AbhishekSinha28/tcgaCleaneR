# Generate PCA functioon

get.pca <- function(data){
  .pca <- function(data, is.log) {
    if(is.log){
      data <- data
    }else{
      data <- log2(data + 1)
    }
    svd <- base::svd(scale(x = t(data), center = TRUE, scale = FALSE))
    percent <- svd$d^2/sum(svd$d^2)*100
    percent <-
      sapply(
        seq_along(percent),
        function(i) {round(percent[i], 1)})
    return(list(
      sing.val = svd,
      variation = percent))
  }
  data.set.names <- names(SummarizedExperiment::assays(data)) #c("HTseq_counts", "HTseq_FPKM", "HTseq_FPKM.UQ")
  pca.cancer.tcga  <- lapply(
    data.set.names,
    function(x){
      .pca(
        data = as.matrix(SummarizedExperiment::assay(data, x)),
        is.log = FALSE)
    })
  names(pca.cancer.tcga) <- data.set.names
  return(pca.cancer.tcga)
}

#df6 <- get.pca(data = df5)
