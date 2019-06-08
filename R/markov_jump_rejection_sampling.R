#' Rejection sample a Markov jump process still active at time = Timepoint
#'
#' @param IntDist Initial distribution of the Markov jump process.
#' @param SubIntMat Sub-intensity matrix of the Markov jump process.
#' @param TimePoint Timepoint at which the Markov jump process must surpass.
#'
#' @return Markov jump process active at time x=Timepoint. Including Sequence of visited states and holding times.
#' @export
#'
#' @examples
RejMJP <- function(IntDist = listAlpha[[6]], SubIntMat = listS[[6]], TimePoint = 0.8 ){
  cnd <-  FALSE
  itr <- 0
  while(!cnd){ # do-while
    sample <- simPH(SubIntMat = SubIntMat, phi0 = IntDist)
    cnd <- sample$tau > TimePoint
    itr <- itr + 1
    # if(itr == 1e4) warning("Iteration surpassed 1e4")
  }
  #Cut-off sample at time x
  RejSample <- list()

  RejSample$t <- c(sample$t[sample$t < TimePoint], TimePoint)

  #seqeunce of visited states, last state prior to abs repeated, for use later e.g. in time spent.
  RejSample$states <- sample$states[sample$t < TimePoint]
  RejSample$states <- append(RejSample$states, tail(RejSample$states,1))

  return(structure(list(t = RejSample$t,
                        states = RejSample$states,
                        IntVec = IntDist,
                        SubIntMat = SubIntMat,
                        TimePoint = TimePoint,
                        itr = itr),
                   class = 'RejMJP')
  )
}
