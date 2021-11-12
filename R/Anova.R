# Anova Function

#' @title Anova test for year and plate effects
#'
#' @description This function is a part of the data analysis functionality of tgcapkg. It helps to perform the Anova test to analyse the variation effect due to time and plate on TCGA Cancer data.
#'
#' @param data S4 data object
#' @param variable character: The predictor variable to \code{lm} model. The variables included are 'Time' and 'Plate'
#' @param is.log logical: Checks if the S4 data has log values. It 'False', it converts data to log scale.
#' @param n.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores.
#'
#' @return A S3 data frame. The output contains the Anova test (F test) scores corresponding to all genes in S4 data object.
#' @export
#'
#' @examples
#' f.test(data = brca.data, variable = "Plate", is.log = FALSE, n.cores = 1)
#' \dontrun{
#'
#' f.test(data = brca.data, variable = "Time", is.log = FALSE, n.cores = 1)
#' }
f.test <- function(data, variable, is.log, n.cores){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  sample.info$ls <- log2(colSums(raw.count))
  raw.count <- as.matrix(raw.count)
  average.exp <- log2(rowMeans(raw.count))
  if(is.log == TRUE){
    raw.count <- raw.count
  }else{
    raw.count <- log2(raw.count + 1)
  }
  if (variable == "Plate"){
    f.test <- parallel::mclapply(
      1:nrow(raw.count),
      function(x) {
        MASS::dropterm(lm(raw.count[x , ] ~ sample.info$PlateId_mda), test = 'F')[c(5:6)]
      }
      , mc.cores = n.cores)
  } else
    if (variable == "Time"){
      f.test <- parallel::mclapply(
        1:nrow(raw.count),
        function(x) {
          MASS::dropterm(lm(raw.count[x , ] ~ sample.info$year_mda), test = 'F')[c(5:6)]
        }
        , mc.cores = n.cores)
    }
  test.values <- data.frame(
    Genes = row.names(raw.count),
    FValue = round(unlist(lapply(f.test, function(x) x$`F Value`[2])), digits = 4) ,
    PValue = unlist(lapply(f.test, function(x) x$`Pr(F)`[2])),
    Adj.PValue = p.adjust(unlist(lapply(f.test, function(x) x$`Pr(F)`[2])), method = 'BH'),
    Mean = round(average.exp, digits = 2)
  )
  return(test.values)
}
