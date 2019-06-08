#'  O-functions
#' Time spent in state `i`
#'
#' @param obj A Markov jump process object with state and jump sequence and intensity matrix.
#' @param State
#'
#' @return
#' @export
#'
#' @examples
#' (a <- RejMJP(TimePoint = 1))
#' OtimeSpent(a, State = 2)
#'
#' (obj <- epsRejMJP(TimePoint = 0.6))
#' OtimeSpent(obj, State = 2)
OtimeSpent <- function(obj, State) {
  if (State %% 1 != 0)
    stop("State must be integer")
  if (!(State %in% 1:dim(obj$SubIntMat)[1]))
    stop("State is not among the transient states.")

  #Calculate time spent in state
  sum(obj$t[which(obj$states == State) + 1] - obj$t[which(obj$states == State)], na.rm = TRUE)
}

#' @describeIn OtimeSpent total number of jump from \eqn{i -> j}
#'
#' @note Note last state in `RejMJP` is not an actual jump
#'
#' @param obj
#' @param i A state in the state space of the Markov jump process.
#' @param j A state in the state space of the Markov jump process.
#'
#' @return The total number of jump from state `i` to `j` of the Markov jump process given by `obj`.
#' @export
#'
#'
#' @examples
#' (obj <- RejMJP(TimePoint = 6))
#' Ojumpsij(obj, i = 2, j = 1)
Ojumpsij <- function(obj, i = 1, j = 2) {
  #Note last state in RejMJP is not an actual jump
  if ((i %% 1) != 0 ||
      (j %% 1) != 0)
    stop("states must be integers")

  numb_states <- dim(obj$SubIntMat)[1]
  if (any(!c(i, j) %in% 1:numb_states))
    stop("i,j are not among states of underlying process")

  len <- length(obj$states)
  SeqStates <- obj$states[-len]

  istates <- which(SeqStates == i)
  sum(SeqStates[istates + 1] == j, na.rm = TRUE)
}

#' @describeIn OtimeSpent total number of jump
#'
#' @note Note last state in `RejMJP` is not an actual jump and first state is not a jump either
#'
#' @param obj
#'
#' @return
#' @export
#'
#'
#' @examples
#' (obj <- RejMJP(TimePoint = 1))
#' OtotalNrJumps(obj = obj)
OtotalNrJumps <- function(obj) {
  #Note last state in RejMJP is not an actual jump and first state is not a jump either
  totalNumbJumps <- length(obj$states) - 2

  totalNumbJumps
}
