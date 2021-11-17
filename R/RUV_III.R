# RUV-III generation

#' @title Generate RUV-III Data Object
#'
#' @description This function is a part of the data analysis functionality of \code{tgcapkg}. It captures both the uses PRPS values from library \code{tgcapkg} combined with row counts and run SVD algorithm (\code{runSVD()}) from \code{BiocSingular} on the combined dataset. The function uses RUV-I algorithm from \code{ruv} as a pre-processing step to RUV-III.
#'
#' @param ruv.data S4 data object for RUV-III: A S4 data object with combined data including the row count from original filtered data using assay \code{HTseq_counts}, prps data for batch and library size. This data needs to be further converted to log scale and transposed.
#' @param ruv.rep S4 data matrix for RUV-III: A S4 data object that has been generated using \code{replicate.matrix} functionality from \code{ruv} package. This helps ruv to identify replicate samples.
#' @param ncg.set logical object: Set of Negative Controlled genes.
#' @param k Integer scalar specifying the number of. Default is NULL. Currently Value 1 represents the library size, 2 represents purity and 3 is time variation.
#' @param eta Gene-wise (as opposed to sample-wise) covariates. A matrix with n columns for \code{ruv::RUV1}. Default is NULL.
#' @param svd_k Integer scalar specifying the number of singular values to return for \code{BiocSingular::runSVD}.Default is 50.
#' @param include.intercept Add an intercept term to eta if it does not include one already for \code{ruv::RUV1}. Default is True.
#' @param BPPARAM A BiocParallelParam object specifying how parallelization should be performed. Default is \code{SerialParam()}
#' @param BSPARAM A BiocSingularParam object specifying the type of algorithm to run. Default is \code{ExactParam()}
#' @param fullalpha To perform RUV-III calculation. Default is NULL.
#' @param return.info logical: Do you want all the information related to RUV-III object. False gives all information whereas True gives only the
#' @param inputcheck logical: Check the inputs to identify if ruv.data contains missing values or infinite values.
#'
#' @return Based on the return.info we get either a S4 list will all information related to RUV-III object or just the RUV-III result.
#' @export
#'
#' @examples
#' \dontrun{
#' get.ruv(ruv.data = ruv.data, ruv.rep = ruv.rep, ncg.set = ncg.set, k=1, BSPARAM = BiocSingular::bsparam(), return.info = TRUE)
#' }

get.ruv <- function(ruv.data, ruv.rep, ncg.set, k = NULL, eta = NULL,
                    svd_k = 50, include.intercept = TRUE,
                    BPPARAM = SerialParam(), BSPARAM = ExactParam(),
                    fullalpha = NULL, return.info = FALSE, inputcheck = TRUE){

  m <- nrow(ruv.data)
  n <- ncol(ruv.data)
  ruv.rep <- ruv::replicate.matrix(ruv.rep)

  tological <- function(ctl, n) {
    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    return(ctl2)
  }

  ctl <- tological(ncg.set, n)

  my_residop <- function(A, B){
    tBB = DelayedArray::t(B) %*% B
    tBB_inv = Matrix::solve(tBB)
    BtBB_inv = B %*% tBB_inv
    tBA = DelayedArray::t(B) %*% A
    BtBB_inv_tBA = BtBB_inv %*% tBA
    return(A - BtBB_inv_tBA)
  }
  ## Check the inputs
  if (inputcheck) {
    if (sum(is.na(ruv.data)) > 0) {
      stop("Input ruv.data contains missing values. This is not supported.")
    }
    if (sum(ruv.data == Inf, na.rm = TRUE) + sum(ruv.data == -Inf, na.rm = TRUE) >
        0) {
      stop("Input ruv.data contains infinities. This is not supported.")
    }
  }
  ## RUV1 is a reprocessing step for RUVIII
  ruv.data <- ruv::RUV1(ruv.data, eta, ctl, include.intercept = include.intercept)
  if (class(BSPARAM) != "ExactParam") {
    svd_k <- min(m - ncol(ruv.rep), sum(ctl), svd_k, na.rm = TRUE)
  } else {
    svd_k <- min(m - ncol(ruv.rep), sum(ctl), na.rm = TRUE)
  }

  ## m represent the number of samples/observations ncol(M)
  ## represent the number of replicates If the replicate matrix
  ## is such that we have more replicates than samples, then
  ## RUV3 is not appropriate, thus, we return the Original input
  ## matrix
  if (ncol(ruv.rep) >= m | k == 0) {
    newY <- ruv.data
    fullalpha <- NULL
  } else {
    if (is.null(fullalpha))
    {
      ## The main RUVIII process Applies the residual operator of a
      ## matrix M to a matrix Y Y0 has the same dimensions as Y,
      ## i.e. m rows (observations) and n columns (genes).
      Y0 <- my_residop(ruv.data, ruv.rep)
      svdObj <- BiocSingular::runSVD(
        x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)

      fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% ruv.data
    }  ## End is.null(fullalpha)
    ###############
    alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]
    W <- ruv.data[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
    newY <- ruv.data - W %*% alpha
  }  ## End else(ncol(M) >= m | k == 0)

  ## If the users want to get all the informations relating to
  ## the RUV, it can be done here.
  if (!return.info) {
    return(newY)
  } else {
    return(list(new.ruv.data = newY, ruv.rep = ruv.rep, fullalpha = fullalpha, W = W))
  }
}
