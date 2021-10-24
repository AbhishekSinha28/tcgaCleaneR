# Anova Function

anova.test <- function(data, variable, is.log, n.cores){
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
  f.test <- data.frame(
    Genes = row.names(raw.count),
    FValue = round(unlist(lapply(f.test, function(x) x$`F Value`[2])), digits = 4) ,
    PValue = unlist(lapply(f.test, function(x) x$`Pr(F)`[2])),
    Adj.PValue = p.adjust(unlist(lapply(f.test, function(x) x$`Pr(F)`[2])), method = 'BH'),
    Mean = round(average.exp, digits = 2)
  )
  return(f.test)
}


#df10 <- anova.test(data = df5, variable = "Plate", is.log = F, n.cores = 5)
#df11 <- anova.test(data = df5, variable = "Time", is.log = F, n.cores = 5)

