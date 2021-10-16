# Vector correlation analysis function

pca.corr <- function(pca.data, data){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  time.years <- fastDummies::dummy_cols(sample.info$year_mda)
  time.years <- time.years[, c(2:ncol(time.years))]
  data.set.names <- names(SummarizedExperiment::assays(data))
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
}
