#  PCs and library size/ purity/ time - Correlation Analysis (regression + vector correlation)

#' @title Correlation Analysis between PCs and variation types
#'
#' @description This function is a part of the data analysis functionality of tgcapkg. It helps user to run regression between bias in TCGA RNA-seq like librarysize and purity with PCs from \code{get.pca}. The output is a linear plot that compares the three \code{assays} in \code{SummarizedExperiment} TGCA Cancer data across n PCs and R-sq. It also runs vector correlation \code{stats::cancor} between time and n PCs with the same linear output explaining variation explained by different variation types.
#'
#' @param pca.data list: PCA output from \code{get.pca}.
#' @param data S4 data object
#' @param type character: The response variable to \code{lm} model. groups included are 'librarysize', 'purity' and 'time'.
#' @param nPCs numeric: Number of PCs that needs to be used for regression
#'
#' @return Linear Plot the compares the correlation between library size (or Purity, time) and PCs across three datasets. When output is stored in a object the user can also access values used to plot the linear graphs.
#' @export
#'
#' @examples
#' \dontrun{
#' pca.corr(pca.data, data = brca.data, type = "purity", nPCs = 10)
#' df <- pca.corr(pca.data, data = brca.data, type = "time", nPCs = 8)
#' df
#' }
pca.corr <- function(pca.data, data, type, nPCs){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library.size <- log2(colSums(raw.count))
  data.set.names <- names(SummarizedExperiment::assays(data))
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  if(type == "librarysize"){
    corr.cancer.tcga <- lapply(
      data.set.names,
      function(x){
        pcs <- pca.data[[x]]$sing.val$u
        tcga.ls.rSquared <- sapply(
          1:nPCs,
          function(y) {
            lm.ls <- summary(lm(library.size ~ pcs[, 1:y]))$r.squared
          })
      })
  } else
    if(type == "purity"){
      corr.cancer.tcga <- lapply(
        data.set.names,
        function(x){
          pcs <- pca.data[[x]]$sing.val$u
          tcga.ls.rSquared <- sapply(
            1:nPCs,
            function(y) {
              purity.ls <- summary(lm(sample.info$purity_HTseq_counts ~ pcs[, 1:y]))$r.squared
            })
        })
    } else
      if(type == "time"){
        time.years <- fastDummies::dummy_cols(sample.info$year_mda)
        time.years <- time.years[, c(2:ncol(time.years))]
        corr.cancer.tcga <-
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
  names(corr.cancer.tcga) <- data.set.names

  # Visualize
  corr.normAssess <- data.frame(
    Raw.counts = corr.cancer.tcga$HTseq_counts,
    FPKM = corr.cancer.tcga$HTseq_FPKM,
    FPKM.UQ = corr.cancer.tcga$HTseq_FPKM.UQ,
    pcs = c(1:nPCs)
  )

  # Visualize
  dataSets.colors <- wesanderson::wes_palette(
    n = 4,
    name = "GrandBudapest1")[c(1,2,4)]
  names(dataSets.colors) <- c(
    'Raw counts',
    'FPKM',
    'FPKM.UQ'
  )

  corr.normAssess.p <- corr.normAssess %>%
    tidyr::pivot_longer(
      -pcs,
      names_to = 'Datasets',
      values_to = 'CorrValue') %>%
    dplyr::mutate(Datasets = replace(
      Datasets,
      Datasets == 'Raw.counts', 'Raw counts')) %>%
    dplyr::mutate(
      Datasets = factor(
        x = Datasets,
        levels = c('Raw counts', 'FPKM', 'FPKM.UQ'))) %>%
    data.frame(.)

  ggplot(corr.normAssess.p, aes(x = pcs, y = CorrValue, group = Datasets)) +
    geom_line(aes(color = Datasets), size = 1) +
    geom_point(aes(color = Datasets), size = 3) +
    xlab('PCs') + ylab (expression("R"^"2")) +
    scale_color_manual(
      values = c(dataSets.colors[1:3]),
      labels = c('Raw counts', 'FPKM','FPKM.UQ')) +
    scale_x_continuous(breaks = (1:nPCs), labels = c('PC1', paste0('PC1:', 2:nPCs)) ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = 'black', size = 1),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 14),
      strip.text.x = element_text(size = 10)
    ) + {if (type == 'purity')
    {
      ggtitle('Purity : PCs Regression Plot')
    } else if (type == 'librarysize'){
      ggtitle('Library Size : PCs Regression Plot')
    } else if (type == 'time'){
      ggtitle('Time : PCs Correlation Plot')
    } }

  #return(corr.normAssess)
}
