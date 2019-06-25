#' Discrete phase-type representation of Urn experiment as given in Example 3.18
#'
#' @param n total number of balls in urn
#' @param d number of black balls in urn
#' @param k k'th black ball
#'
#' @return Phase type representation (`pi_kn` , `T_kn`) as a list.
#' @export
#'
#' @examples
#'
#' Urn3(100, 25, 10)
Urn3 <- function(n, d, k){
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()

  #* Add an error for n , number of balls
  if ( n%%1 != 0 | n <= 0)
    ArgumentCheck::addError(
      msg = "n must be positive integers",
      argcheck = Check
    )

  #* Add an error for d, number of black balls
  if ( d%%1 != 0 | d > n | d <= 0)
    ArgumentCheck::addError(
      msg = "d must be integers, 0 < d <= n",
      argcheck = Check
    )

  #* Add an error for k, the draw of black ball
  if ( k%%1 != 0 | k > d | k <= 0)
    ArgumentCheck::addError(
      msg = "k must be integers, 0 =< k <= d",
      argcheck = Check
    )


  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  ############################### Checks Done

  n <- as.integer(n)
  d <- as.integer(d)
  k <- as.integer(k)


  #Initial vector
  #start in state 1 , out of total state space
  pi_kn <- rep(0, (n-d+1)*k) #Initial of multidim n-process
  pi_kn[1] <- 1

  #Special Case: n = d, all balls are black.
  #Simply return k-matrix with 1's on superdiagonal.
  if(n == d){
    T_kn <- matrix(0, nrow = k, ncol = k)
    #fill super-diagonal with 1s
    if(k >= 2){
      for(i in seq_along(k-1))
        T_kn[i,i+1] <- 1
    }

    return(list(T_kn = T_kn,
                pi_kn = pi_kn))
  }


  #T_(j:n) subtransition matrix
  #T_jn function to create Blocks
  T_jn <- function(n, d, j=1){
    #Matrix of dim n-d+1
    out <- matrix(0, ncol = n-d+1, nrow = n-d+1)
    numerator <- (n-(j-1)-(d-j+1)):1
    denominator <- (n-(j-1)):(d-(j-1)+1)

    for(i in 1:(dim(out)[1]-1)){
      out[i,i+1] <- numerator[i]/denominator[i]
    }
    out
  }
  #T0_jn function to create the exit matrices
  T0_jn <- function(n, d , j = 1){
    #matrix of dim n-d+1
    out <- matrix(0, ncol = n-d+1, nrow = n-d+1)
    numerator <- d-(j-1)
    denominator <- (n-(j-1)):(d-(j-1))

    diag(out) <- numerator/denominator
    out
  }


  #Create zero-filled T_k:n matrix and fill out
  #Fill in the Diag Blocks
  listT_jn <- lapply(1:k, function(x) T_jn(n = n , d = d, j = x))
  T_kn <- Reduce(listT_jn, f = function(x,y) Matrix::bdiag(x,y))
  T_kn <- as.matrix(T_kn)


  #Fill in S^0_kn blocks
  if(k >= 2){
    dimBlock <- dim(listT_jn[[1]])[1]
    CumDim <- cumsum(unlist(lapply(listT_jn, function(x)dim(x)[1])))

    for(m in 1:(k-1)){
      T0jn <- T0_jn(n = n,d = d, j = m)

      T_kn[ (CumDim[m] - dimBlock + 1):CumDim[m], (CumDim[m] + 1):(CumDim[m]+dimBlock)] <- T0jn

    }
  }

  return(list(T_kn = T_kn,
              pi_kn = pi_kn))

}
