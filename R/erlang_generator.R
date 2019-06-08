#' Erlang Feedback example
#' General Erlang (feedback) generator
#'
#' @param k number of phases in the underlying process.
#' @param rate rate parameter. Can also be a k vector.
#' @param feedback Feedback probability parameter, default value is zero.
#'
#' @return
#' @export
#'
#' @examples
erlang_ph <- function(k = 10, rate = 1:10, feedback = NULL){
  if( length(rate) != 1 && length(rate) != k   ) stop("Length of rate must equal either one or number of phases")
  if( any(rate <= 0) ) stop("Rates must all be positive")

  Tmat <- matrix(0,k,k)

  diag(Tmat) <- -rate
  if(length(rate) ==1 ){
    diag(Tmat[1:(k-1),2:k]) <- rate
  } else {
    diag(Tmat[1:(k-1),2:k]) <- rate[1:(k-1)]
  }

  if(!is.null(feedback) ){
    if( feedback < 0 || feedback > 1) stop("feedback prop must be between 0 and 1")
    Tmat[k,1] <- feedback*(-Tmat[k,k])
  }
  InitVec <- c(1,rep(0,k-1))
  return(list(Tmat=Tmat,
              InitVec = InitVec))
}
