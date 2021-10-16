#  PCs and library size/ purity - linear regression

pca.reg <- function(pca.data, data, type){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library.size <- log2(colSums(raw.count))
  data.set.names <- names(SummarizedExperiment::assays(data))
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  nPCs <- 10
  if(type == "librarysize"){
    lreg.cancer.tcga <- lapply(
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
      lreg.cancer.tcga <- lapply(
        data.set.names,
        function(x){
          pcs <- pca.data[[x]]$sing.val$u
          tcga.ls.rSquared <- sapply(
            1:nPCs,
            function(y) {
              purity.ls <- summary(lm(sample.info$purity_HTseq_FPKM ~ pcs[, 1:y]))$r.squared
            })
        })
    }
  names(lreg.cancer.tcga) <- data.set.names

  # Visualize
  lnreg.normAssess <- data.frame(
    Raw.counts = lreg.cancer.tcga$HTseq_counts,
    FPKM = lreg.cancer.tcga$HTseq_FPKM,
    FPKM.UQ = lreg.cancer.tcga$HTseq_FPKM.UQ,
    pcs = c(1:10)
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

  lnreg.normAssess <- lnreg.normAssess %>%
    tidyr::pivot_longer(
      -pcs,
      names_to = 'Datasets',
      values_to = 'ls') %>%
    dplyr::mutate(Datasets = replace(
      Datasets,
      Datasets == 'Raw.counts', 'Raw counts')) %>%
    dplyr::mutate(
      Datasets = factor(
        x = Datasets,
        levels = c('Raw counts', 'FPKM', 'FPKM.UQ'))) %>%
    data.frame(.)

  ggplot(lnreg.normAssess, aes(x = pcs, y = ls, group = Datasets)) +
    geom_line(aes(color = Datasets), size = 1) +
    geom_point(aes(color = Datasets), size = 3) +
    xlab('PCs') + ylab (expression("R"^"2")) +
    scale_color_manual(
      values = c(dataSets.colors[1:3]),
      labels = c('Raw counts', 'FPKM','FPKM.UQ')) +
    scale_x_continuous(breaks = (1:10), labels = c('PC1', paste0('PC1:', 2:10)) ) +
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
    )
}

#pca.reg(pca.data = df5, data = df4, type = "librarysize")
