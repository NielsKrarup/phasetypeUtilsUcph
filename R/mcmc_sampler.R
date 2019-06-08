#' Markov Chain Monte Carlo Estimator of functional of underlying Markov jump process underlying a Phase-type distribution conditioned on the absorption time being equal to x = TimePoint.
#'
#' @param pi The initial distribution of the underlying Markov jump process.
#' @param SubIntMat The sub-intensity matrix of the underlying Markov jump process.
#' @param N Number of steps in the MCMC algorithm.
#' @param BurnIn Burn In.
#' @param TimePoint Time of absorption conditioned on.
#' @param OFun The functional of the underlying process.
#' @param ...
#'
#' @return Estimated mean of function `Ofun` applied to the underlying process.
#' @export
#'
#' @examples
#' MCMC_PH(pi = listAlpha[[4]], SubInt = listS[[4]], N = 1e4, BurnIn = 1e2, TimePoint = 0.8, OFun = OtimeSpent, State = 5)
MCMC_PH <-
  function(pi = ErInt,
           SubIntMat = ErMat,
           N = 1e2,
           BurnIn = NULL,
           TimePoint = 0.8,
           OFun = OtimeSpent,
           ...) {
    #The general state Markov chain.
    #FIXME: use `vector(mode="list", N)` as initializer.
    Zn <- list()
    Zn[[1]] <-
      RejMJP(IntDist = pi,
             SubIntMat = SubIntMat,
             TimePoint = TimePoint)

    #running the MH

    #Draw U beforehand for speed
    Uvec <- runif(n = N)

    for (i in 1:(N - 1)) {
      #Proposal: process surpassing Timepoint
      Y <-
        RejMJP(IntDist = pi,
               SubIntMat = SubIntMat,
               TimePoint = TimePoint)

      ratio <-
        -rowSums(Y$SubIntMat)[tail(Y$states, 1)] / -rowSums(Zn[[i]]$SubIntMat)[tail(Zn[[i]]$states, 1)]
      if (is.nan(ratio)) {
        ratio <- 0
      }

      if (Uvec[i] < min(1, ratio)) {
        Zn[[i + 1]] <- Y
      } else{
        Zn[[i + 1]] <- Zn[[i]]
      }
    }

    #BurnIn
    if (!is.null(BurnIn)) {
      Zn <- Zn[-c(1:BurnIn)]
    }

    #Calculate functionals
    # Sampler, untill surpassed x = Timepoint

    #Importance sampler mean 1/N sum_1^N Weights*O(n)

    sample <- unlist(lapply(
      X = Zn,
      FUN = function(x)
        OFun(obj = x, ...)
    ))


    sample_mean <- mean(sample)
    sample_var  <- var(sample)
    sample_sd <- sd(sample)


    return(list(
      sample_mean = sample_mean,
      sample_var = sample_var,
      sample_sd = sample_sd
    ))

  }
