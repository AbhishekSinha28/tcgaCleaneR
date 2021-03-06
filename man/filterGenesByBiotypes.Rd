% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_genes.R
\name{filterGenesByBiotypes}
\alias{filterGenesByBiotypes}
\title{Gene Filter}
\usage{
filterGenesByBiotypes(data,gene.type)
}
\arguments{
\item{data}{SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.}

\item{gene.type}{A character vector of items.}
}
\value{
S4 data object
}
\description{
This function is a part of the data wrangling functionality of \code{tcgaCleaneR}.
It allows user to input the \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) and the required
genes to filter data based on genes.
}
\examples{
filterGenesByBiotypes(data=brca.data,gene.type=c("protein.coding"))
\dontrun{
filterGenesByBiotypes(data=brca.data,gene.type=c("protein.coding","snRNA"))
}
}
