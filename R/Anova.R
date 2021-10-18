# Anova Function

.Ftest <- function(data, variable, is.log, n.cores){
  average.exp <- log2(rowMeans(data))
  if(is.log == TRUE){
    data <- data
  }else{
    data <- log2(data + 1)
  }
  f.test <- parallel::mclapply(
    1:nrow(data),
    function(x) {
      MASS::dropterm(lm(data[x , ] ~ variable), test = 'F')[c(5:6)]
    }
    , mc.cores = n.cores)
  f.test <- data.frame(
    Genes = row.names(data),
    FValue = round(unlist(lapply(f.test, function(x) x$`F Value`[2])), digits = 4) ,
    PValue = unlist(lapply(f.test, function(x) x$`Pr(F)`[2])),
    Adj.PValue = p.adjust(unlist(lapply(f.test, function(x) x$`Pr(F)`[2])), method = 'BH'),
    Mean = round(average.exp, digits = 2)
  )
  return(f.test)
}
