#' Importance Sampler of process underlying a Phase-type distribution conditioned on absorption time.
#'
#' @param pi Initial distribution of underlying process.
#' @param SubIntMat Sub-intensity matrix of underlying process.
#' @param N Number of total iterations.
#' @param TimePoint Time of absorption to be conditioned on.
#' @param OFun Functional applied to underlying process.
#' @param ... Additional argument passed on to `OFun`.
#'
#' @return
#' @export
#'
#' @examples
#' ImpSampler(State = 1)
ImpSampler <-
  function(pi = listAlpha[["pi4"]],
           SubIntMat = listS[["S4"]],
           N = 1e2,
           TimePoint = 0.8,
           OFun = OtimeSpent,
           ...) {
    #Should be SubIntMat Mat
    if (any(diag(SubIntMat) == 0))
      stop("must be SubIntMat Matrix supplied.")
    if (length(pi) != dim(SubIntMat)[1])
      stop("Check dimension of IntDist and/or SubIntMat")
    if (sum(pi) != 1)
      stop("int distribution must sum to one")


    # Sampler, untill surpassed x = Timepoint
    tmplist <-
      replicate(n = N, expr = list(
        RejMJP(
          IntDist = pi,
          SubIntMat = SubIntMat,
          TimePoint = TimePoint
        )
      ))

    #Importance sampler mean 1/N sum_1^N Weights*O(n)
    sample <- unlist(lapply(
      X = tmplist,
      FUN = function(x) {
        -rowSums(x$SubIntMat)[tail(x$states, 1)] * OFun(obj = x, ...)
      }
    )) * rowSums(pi %*% Matrix::expm(SubIntMat * TimePoint)) /
      as.numeric(pi %*% Matrix::expm(SubIntMat *
                               TimePoint) %*% (-rowSums(SubIntMat)))


    return(list(
      N = N,
      sample_mean = mean(sample),
      sample_var = var(sample)
    ))
  }

