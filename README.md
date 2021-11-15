
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tgcapkg

<!-- badges: start -->

![R-CMD-check](https://github.com/AbhishekSinha28/tgcapkg/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

The goal of `tgcapkg` is to provide a user-friendly R package to help
Bioinformaticians with easy access to a tool that can perform Data
Wrangling and Data Analysis on TCGA Pan Cancer Dataset. The package
contains a subset of the original TCGA Breast Cancer Data collected from
[TCGA](https://www.cbioportal.org/datasets). The package contains a
detailed set of functionalities allowing user to identify and handle
unwanted variations in the TCGA datasets.

## Installation

You can install the development version of tgcapkg from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AbhishekSinha28/tgcapkg", ref="master", auth_token = "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")
```

## TCGA Functionality

This is a quick walk trough the `tgcapkg` functionalities. For the
detailed information on the Data Wrangling and Data Analysis functions
and arguments of `tgcapkg` package, files can be refereed. Also you can
consider looking at the vignette to get a detailed walk through about
the package functionalties.

``` r
library(tgcapkg)
```

# Data

``` r
data("brca.data")
```

``` r
brca.data
#> class: SummarizedExperiment 
#> dim: 100 1222 
#> metadata(0):
#> assays(3): HTseq_counts HTseq_FPKM HTseq_FPKM.UQ
#> rownames(100): ENSG00000000003 ENSG00000000005 ... ENSG00000005469
#>   ENSG00000005471
#> rowData names(40): gene_id.v gene_id ... keep.cancer keep.normal
#> colnames(1222): TCGA-A8-A06Z-01A-11R-A00Z-07
#>   TCGA-A8-A08F-01A-11R-A00Z-07 ... TCGA-BH-A0H9-11A-22R-A466-07
#>   TCGA-E9-A54Y-01A-11R-A466-07
#> colData names(4115): Sample sample.id.b_mda ... purity_HTseq_FPKM
#>   purity_HTseq_FPKM.UQ
```

# Data Wrangling

## Gene Filter

``` r
filtered.data <- gene.filter(data=brca.data,gene.type=c("protein_coding"))
```

## Removing lowly expressed genes

``` r
filtered.data1 <- low.genes.filter(data=filtered.data,gene_count = 20,sample_size = 200)
```

## Purity Filter - Filter Samples based on Tumor Purity

``` r
filtered.data2 <- purity.filter(data= filtered.data1,purity_cutoff= 0.496)
```

## Library Size Filter

### Determine Library Size

``` r
library.size(data = filtered.data2, plot_type = "Scatterplot")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" style="display: block; margin: auto;" />

### Filter samples based on library size

``` r
filtered.data3 <- library.size.filter(data = filtered.data2, ls_cutoff = 17.5)
```

# Data Analysis

## Study Design Plot

The idea behind Study Design plot is to present the summarized
information about the filtered dataset using HeatMaps.

``` r
study.design(data = filtered.data3)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" style="display: block; margin: auto;" />

## PCA

### Generate PCA

The principal components (in this context also called singular vectors)
of the sample × transcript array of log-counts are the linear
combinations of the transcript measurements having the largest, second
largest, third largest, etc. variation, standardized to be of unit
length and orthogonal to the preceding components. Each will give a
single value for each sample.

``` r
# Is data input for PCA logical
is.logical(filtered.data3)
#> [1] FALSE
```

``` r
pca_data <- get.pca(data = filtered.data3, nPcs = 7, is.log = FALSE)
```

``` r
# Generated 7 PCs for first two genes
head(pca_data[['HTseq_counts']]$sing.val$u,2)
#>                                      [,1]        [,2]        [,3]        [,4]
#> TCGA-A8-A06Z-01A-11R-A00Z-07 -0.008247742 -0.03800256 0.074868654  0.08068739
#> TCGA-AN-A03Y-01A-21R-A00Z-07  0.001517381 -0.01297216 0.008098276 -0.01609687
#>                                      [,5]       [,6]        [,7]
#> TCGA-A8-A06Z-01A-11R-A00Z-07 -0.017982618 0.03732317  0.02341218
#> TCGA-AN-A03Y-01A-21R-A00Z-07 -0.001079871 0.04753709 -0.04164148
```

### Plot PCA

Once we have the PCs generated using the PCA function the next step is
to visualise those PCs with respect to the sample features like Time,
Tissue, Plate etc., to identify any unwanted variation by identifying
patterns in the plots by feature.

``` r
library(ggplot2)
library(cowplot)
pca.plot.data <- pca.plot(pca.data = pca_data, data = filtered.data3, group = "Time", plot_type = "DensityPlot", npcs = 3)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-14-2.png" width="100%" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-14-3.png" width="100%" style="display: block; margin: auto;" />

### PCs correlation with unwanted variations

``` r
library(tidyverse)
corr_data <- pca.corr(pca.data = pca_data, data = filtered.data3, type = "purity", nPCs = 7)
corr_data
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" style="display: block; margin: auto;" />
### Gene correlation

``` r
gene_corr_data <- gene.corr(data = filtered.data3, is.log = FALSE, type = "librarysize", cor.method = 'spearman', n.cores = 2)
#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties

#> Warning in cor.test.default(x = expr.data[x, ], y = variable, method = method):
#> Cannot compute exact p-value with ties
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" style="display: block; margin: auto;" />

``` r
head(gene_corr_data)
#>             genes      rho       pvalue   adj.pvalue
#> 1 ENSG00000000003 0.070094 1.462921e-01 1.614258e-01
#> 2 ENSG00000000005 0.085557 7.601159e-02 8.687039e-02
#> 3 ENSG00000000419 0.497268 2.656644e-28 1.275189e-27
#> 4 ENSG00000000457 0.455709 1.730827e-23 6.923306e-23
#> 5 ENSG00000000460 0.456292 1.496799e-23 6.247507e-23
#> 6 ENSG00000000938 0.213431 7.833556e-06 1.212938e-05
```

# Anova Test

``` r
anova.data <- anova_test(data = filtered.data3, variable = "Time", is.log = FALSE, n.cores = 2)
head(anova.data)
#>                           Genes  FValue       PValue   Adj.PValue  Mean
#> ENSG00000000003 ENSG00000000003  1.9074 1.082485e-01 1.404305e-01 11.80
#> ENSG00000000005 ENSG00000000005  1.5990 1.735952e-01 2.057425e-01  5.19
#> ENSG00000000419 ENSG00000000419 11.0646 1.512727e-08 6.915321e-08 11.35
#> ENSG00000000457 ENSG00000000457  9.4292 2.602528e-07 9.609334e-07 11.09
#> ENSG00000000460 ENSG00000000460  6.5534 3.993934e-05 1.161872e-04  9.92
#> ENSG00000000938 ENSG00000000938  4.5189 1.384707e-03 3.021178e-03  8.21
```
