#' Columns of matrix as list of columns
#'
#' @param x
#'
#' @return
#'
#' @examples
cols_to_list <-
  function(x) {
    lapply(seq_len(ncol(x)), function(i)
      x[, i])
  }

# cols_to_list(matrix(1:9, ncol = 3))

#' The dimension of the k'th order stat of n PH
#'
#' @param k the order.
#' @param n total number of phase-types from which to calculate dimension of kth order statistic.
#' @param pvec vector of dimensions of n phase-types. If all equal can be single number.
#'
#' @return Dimension of the phase-type representation of the kth order statistic of n.
#' @export
#'
#' @examples
dimOrderStat <- function(k, n, pvec) {
  if (k > n | k < 1)
    stop('k must be between 1 and n')
  if (length(pvec) != n) {
    pvec <- rep(pvec, n)
  }

  sum(sapply(
    1:k,
    FUN = function(x)
      sum(apply(combn(pvec, n - x + 1), 2, prod))
  ))
}


#' Kronecker SUM
#'
#' @param A Square matrix.
#' @param B Square matrix.
#'
#' @return The Kronecker Sum between matrices A and B.
#' @export
#'
#' @examples sumKronecker(listS[[1]], listS[[2]])
sumKronecker <- function(A, B) {
  if (dim(A)[1] != dim(A)[2] |
      dim(B)[1] != dim(B)[2])
    stop('Must be square matrices')
  b <- dim(B)[1]
  a <- dim(A)[1]

  Ib <- diag(1, b)
  Ia <- diag(1, a)

  kronecker(A, Ib) + kronecker(Ia, B)
}
