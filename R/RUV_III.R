# RUV-III generation

get.ruv <- function(ruv.data, ruv.rep, ncg.set, k = NULL, eta = NULL,
                    svd_k = 50, include.intercept = TRUE, average = FALSE,
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
