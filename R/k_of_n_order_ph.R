#' Generate the sub-transition matrix of the k:n order statistic given \eqn{X_1,..., X_n} phase-types.
#'
#' @param listAlpha list of initial distributions of the n phase-type distributions.
#' @param listS  list of sub-intensity matrices of the n phase-type distributions.
#' @param k the order statistic of \eqn{n}.
#'
#' @return Sub-intensity matrix and initial vector of phase-type representation of kth order statistic. Also dimensions of each block matrix.
#' @export
#'
#' @examples
knOrderPH <- function(listAlpha, listS, k) {
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()

  #* Add an error if lists not equal length
  if (length(listAlpha) != length(listS))
    ArgumentCheck::addError(msg = "Number of sub-matrices and initial vectors do not match. They should be equal",
                            argcheck = Check)

  #* Add an error if type is not lists
  if (typeof(listS) != 'list' | typeof(listAlpha) != 'list')
    ArgumentCheck::addError(msg = "Input should be of type 'lists'",
                            argcheck = Check)

  n <- length(listS)

  #* Add an error if not 1 <= k <= n or k is not integer
  if (k %% 1 != 0 | k < 1 | k > n)
    ArgumentCheck::addError(msg = "k must be integer between 1 and n",
                            argcheck = Check)

  k <- as.integer(k)

  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  if (n == 1) {
    return(list(T_kn = listS[[1]],
                pi_kn = listAlpha[[1]]))
  }


  #Vector of dimensions of each PH
  pvec <- sapply(listAlpha, length)
  #Dimension of knOrder
  dimkn <-
    dimOrderStat(k = k, n = n, pvec = pvec) #dimension of k'th order stat PH

  pi_kn <-
    Reduce(listAlpha, f = kronecker) #Initial of multidim n-process

  if (dimkn != length(pi_kn)) {
    NumAppendZero <- dimkn - length(pi_kn)
    pi_kn <- c(pi_kn, rep(0, NumAppendZero))
  }


  #T_kn subtransition matrix

  #S_kn function to create Blocks of Diag of Kronecker Sums of subset, S_k:n in book ####
  S_kn <- function(listSubMat = listS, k) {
    #Create list of kronecker sums

    if (k == 1) {
      #One block of Kronecker sum of all S1,...,Sn
      tmp1ListKroSums <-
        list(Reduce(listSubMat[combn(1:n, n - k + 1)], f = sumKronecker))
    } else if (k == n) {
      #n Blocks diag of S1,S2,...,Sn
      tmp1ListKroSums <- listSubMat
    } else
      tmp1ListKroSums <-
        apply(
          X = as.matrix(combn(1:n, n - k + 1)),
          MARGIN = 2,
          FUN = function(x)
            Reduce(listSubMat[x], f = sumKronecker)
        )

    #diag block matrix of kronecker SUMS
    s_kn <- Reduce(
      tmp1ListKroSums,
      f = function(x, y)
        Matrix::bdiag(x, y)
    )
    as.matrix(s_kn)
  }

  #Fill out Diag Block Part of T_kn

  #EXIT MATRICES: S^0_kn
  #I matrices, same dimension as corresponding exit vectors
  listI <- lapply(listS, function(x)
    diag(1, dim(x)))
  listsexit <- lapply(listS, function(x)
    as.matrix(-rowSums(x)))
  listIsexit <- c(listI, listsexit)

  #Helper Fun in S0_kn ####
  entrySzero <- function(aindex, bindex) {
    aindex <- as.vector(aindex)
    bindex <- as.vector(bindex)
    if (all(c(bindex) %in% aindex)) {
      sexitIndex <- setdiff(aindex, bindex)

      IsexitIndex <- aindex
      IsexitIndex[IsexitIndex == sexitIndex] <-
        n + sexitIndex #To select the s-exit vector

      out <- Reduce(listIsexit[IsexitIndex], f = kronecker)
    }
    else {
      atmp <- lapply(listS[c(aindex)], dim)
      adimtmp <- Reduce(
        atmp,
        f = function(x, y)
          x * y
      )

      btmp <- lapply(listS[c(bindex)], dim)
      bdimtmp <- Reduce(
        btmp,
        f = function(x, y)
          x * y
      )
      out <- matrix(0, nrow = adimtmp[1], ncol = bdimtmp[1])
    }
    out
  }

  #Function to calculate the Exitmatrices of S^0_k:n ####
  S0_kn <- function(listSubMat, k) {
    fromStateSpace <- as.matrix(combn(1:n, n - k + 1))
    toStateSpace <- as.matrix(combn(1:n, n - k))

    S0_knExitBlock <- NULL
    for (i in 1:ncol(fromStateSpace)) {
      #Using helper fun to get each exitblock, then cbind them.
      ExitBlocks <-
        lapply(cols_to_list(toStateSpace) , function(b)
          entrySzero(aindex = fromStateSpace[, i] , bindex = b))
      tmp <- Reduce(f = cbind, ExitBlocks)
      S0_knExitBlock <- rbind(S0_knExitBlock, tmp)
    }
    S0_knExitBlock
  }


  #Create zero-filled T_k:n matrix and fill out
  #Fill in the Diag Blocks, made with S_kn()
  listS_kn <-
    lapply(1:k, function(x)
      S_kn(listSubMat = listS, k = x))

  T_kn <- Reduce(
    listS_kn,
    f = function(x, y)
      Matrix::bdiag(x, y)
  )
  T_kn <- as.matrix(T_kn)
  dim(T_kn)

  dimSknBlocks  <-
    lapply(
      listS_kn,
      FUN = function(x)
        dim(x)[1]
    ) #Dimension of used S_kn, square blocks
  cumDim <- cumsum(unlist(dimSknBlocks))

  #Fill in S^0_kn blocks
  if (k > 1) {
    for (j in 1:(k - 1)) {
      S0jnMat <- S0_kn(listSubMat = listS, k = j)
      dimS0jn <-
        dim(S0jnMat) # Dimension of S^0_j:n matrix, not necessarily square!

      T_kn[(cumDim[j] - dimS0jn[1] + 1):cumDim[j][1],
           (cumDim[j] + 1):(cumDim[j] +
                              dimS0jn[2])] <- S0jnMat
    }
  }

  return(list(
    T_kn = T_kn,
    pi_kn = pi_kn,
    cumDimIndex = cumDim
  ))

}
