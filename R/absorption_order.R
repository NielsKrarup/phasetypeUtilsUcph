#' Quantities used in calculation of absorption order probabilities.
#'
#' @param listS list of sub-intensity matrices of the Phase-types.
#' @param listAlpha list of initial distribution vectors of the Phase-types.
#'
#' @return The vector `pi, T, D` of theorem 4.23.
#' @export
#'
#' @examples
AbsOrderPH <- function(listS, listAlpha){

  if(length(listS) != length(listAlpha)) stop("Must supply lists of Initial vectors and Sub intensities matrices, of same length")
  if( any(unlist(lapply(listS, function(x)dim(x)[1])) != unlist(lapply(listAlpha, length))))
    stop("The dimension of initial vector and sub intensity matrix should be equal")
  if(length(listS) < 2) stop("n should be at least 2")

  #####Helper lists
  n <- length(listS)
  #List of dimensions
  listDim <- lapply(listS, function(x) dim(x)[1] )
  #I matrices, same dimension as corresponding exit vectors
  listI <- lapply(listS, function(x) diag(1, dim(x) ) )
  #list of corresponding exit vectors
  listsexit <- lapply(listS, function(x) as.matrix(-rowSums(x)))
  #Combined list of I and s
  listIsexit <- c(listI, listsexit)


  #### Helper Funs ####
  #Helper Fun to get dimension of lower order state space, i.e. for Index = 1,2,4 -> dim of (1,2) + 1 + 2 , (1,4) + 1 + 4, (2,4) + 2+ 4
  Dimension_lower_order <- function(Index = c(1,2,3)){

    n <- length(Index)
    #Dimension of kronecker sum after removing k
    if(length(Index) == 2) {
      sum(unlist(listDim[Index]))
    }
    else {
      #Kronecker sums
      sum( apply(as.matrix(combn(x = Index, n-1)), MARGIN = 2, FUN = function(x) prod(unlist(listDim[x])) ) ) +
        sum(apply(as.matrix(combn(x = Index, n-1)), MARGIN = 2, FUN = function(x) Dimension_lower_order(Index =  x) ))
    }

  }

  #### Helper Fun . Given list of current indices and the index to go, creates the exit matrix. ####
  FirstExitBlock <- function(Index, outIndex){
    if(  !(outIndex %in% Index) ) stop("OutIndex should be in Index")

    Index[Index == outIndex] <- length(listS) + outIndex #To select the s-exit vector

    out <- Reduce(listIsexit[Index], f = kronecker)

    out
  }

  #Since we should only go one layer down, the tail of exit matrix should be zero filled
  SecondZeroExitBlock <- function(Index, outIndex){
    if( !(outIndex %in% Index)) stop("OutINdex must be one of the index")
    if(length(Index) == 2) return(NULL)
    nrows <- prod(unlist(listDim[Index]))
    ncols <- Dimension_lower_order(Index = setdiff(Index,outIndex))
    out <- matrix(0, nrow = nrows, ncol = ncols)
    out
  }


  #Whole block of exit intensities
  ExitBlock <- function(Index = c(1,2,3)){
    if(length(Index) == 1) return( listsexit[[Index]])

    #function to concatenate matrices together
    addBlocks <- function(Index, outIndex){

      FirstExitBlock <- FirstExitBlock(Index = Index, outIndex = outIndex)
      SecondZeroExitBlock <- SecondZeroExitBlock(Index = Index, outIndex = outIndex)

      if(!is.null(SecondZeroExitBlock)){
        if(nrow(FirstExitBlock)!=nrow(SecondZeroExitBlock)) stop("Not same row dimension")
      }

      cbind(FirstExitBlock,SecondZeroExitBlock)
    }

    #Note due to lexicographical ordering, the last Index is removed first.
    listExitBlocks <- lapply(X = rev(Index), FUN = function(x) addBlocks(Index = Index, outIndex = x))
    Reduce(f = cbind, x = listExitBlocks)
  }

  #Creates the diag block Matrix
  BlockGenerator <- function(MatListIndex = c(1,2,3), listS, listAlpha){
    #MatListIndex: the indices of the subintensity matrices in question
    #If only two indices then return subint mat for max.
    n <- length(MatListIndex)

    if(length(MatListIndex)  == 2){
      #Simply Max order stat of n = 2 Phase-types S1, S2
      knOrderPH(listAlpha = listAlpha[MatListIndex], listS = listS[MatListIndex], k = 2)$T_kn
    } else {
      Matrix::bdiag(Reduce(f = sumKronecker, x = listS[MatListIndex]),
            Reduce(f = Matrix::bdiag, x = apply(X = as.matrix(combn(x = MatListIndex, m = n-1)), MARGIN = 2,
                                        FUN = function(b) BlockGenerator(MatListIndex = b, listS = listS,listAlpha = listAlpha))) )
    }
  }



  #### Inserting Exit Matrices.####
  # Find all, ordered combinations exit matrices for appending. given n variables.
  AllExitBlocksIndex <- function(OrderLayer = 3){
    n_counter <- length(listS)

    if(length(listS) == OrderLayer) return(matrix(1:n,ncol = 1))

    #layer n-1
    A <- combn(1:length(listS),  n_counter-1)
    #carefull, counter is only for going down in number of variables
    n_counter <- n_counter-1

    while(n_counter > OrderLayer){
      A <- lapply( cols_to_list(A),FUN = function(x) combn(x,n_counter-1) )
      A <- Reduce(f = cbind, x = A)
      n_counter <- n_counter-1
    }
    A <- as.matrix(A)
    A
  }

  #OrderLayer is from TOP i.e. n first
  Dim_AllExitBlocks <- function(OrderLayer = 3){
    n_counter <- length(listS)

    if(length(listS) == OrderLayer){
      return(as.matrix(dim(ExitBlock(Index = 1:length(listS)))))
    }

    #layer n-1
    A <- combn(1:length(listS),  n_counter-1)
    #carefull, counter is only for going down in number of variables
    n_counter <- n_counter-1

    while(n_counter > OrderLayer){
      A <- lapply( cols_to_list(A),FUN = function(x) combn(x,n_counter-1) )
      A <- Reduce(f = cbind, x = A)
      n_counter <- n_counter-1
    }

    Dim_AllExitBlocks <-   apply(A, 2, FUN = function(x) dim(ExitBlock(Index = as.vector(x) )))
    Dim_AllExitBlocks
  }



  #Get list of All exit-Index for Exit Matrices
  list_IndexAllExitBlocks <- lapply(X = length(listS):1, FUN = function(b) AllExitBlocksIndex(OrderLayer = b) )

  #Get list of ALL Dimensions of Exit matrices in order
  list_DimAllExitBlocks <- lapply(X = length(listS):1, FUN = function(b) Dim_AllExitBlocks(OrderLayer = b))

  #For the rest of the entries a function will insert the UpperLeftRowCol matrix
  UpdaterFun <- function(list_DimAllExitBlocks = list_DimAllExitBlocks, layer_from_top = 3){

    list_Number_of_ExitMatrices <- lapply(list_DimAllExitBlocks, ncol)
    Mat <- matrix(data = NA, nrow = 2, ncol = list_Number_of_ExitMatrices[[layer_from_top]])

    if(list_Number_of_ExitMatrices[[layer_from_top-1]]*(n:1)[layer_from_top-1] !=  list_Number_of_ExitMatrices[layer_from_top])
      stop('Groupsize * number of groups != total number')

    #At layer_from_top the group size is n:3[layer_from_top-1]
    #In layer: layer_from_top
    #loop over each block
    if(layer_from_top != length(listS)){
      number_groups <-list_Number_of_ExitMatrices[[layer_from_top-1]]
      for(k in 1:number_groups){
        #inside each block, number_within_each_group
        number_within_each_group <- (n:1)[layer_from_top-1]
        for(i in 1:number_within_each_group){
          #row, each group should start from corresponding prior k'th row,col
          #Note 1st in each group k should add to prior k'th ExitBlock (row,col) nrow of k'th Exitbloxks
          #Else the nrow + ncol of the prior exit block should be added to row of i
          Mat[1,(k-1)*number_within_each_group + i] <- if(i == 1){
            List_Mat_UpperLeftRowCol[[ layer_from_top -1 ]][1,k] + list_DimAllExitBlocks[[ layer_from_top -1]][1,k]
          }else Mat[1,(k-1)*number_within_each_group + i -1] + sum(list_DimAllExitBlocks[[ layer_from_top]][ , (k-1)*number_within_each_group + i -1])
          #col
          Mat[2,(k-1)*number_within_each_group + i] <- if(i == 1){
            List_Mat_UpperLeftRowCol[[ layer_from_top -1 ]][2,k] + list_DimAllExitBlocks[[ layer_from_top ]][1,(k-1)*number_within_each_group + i]
          } else Mat[2,(k-1)*number_within_each_group + i -1] + list_DimAllExitBlocks[[ layer_from_top]][2,(k-1)*number_within_each_group + i -1] + list_DimAllExitBlocks[[ layer_from_top]][1,(k-1)*number_within_each_group + i]

        }
      }
    }

    #Last layer special case
    if(layer_from_top == length(listS)){
      number_groups <-list_Number_of_ExitMatrices[[layer_from_top-1]]
      for(k in 1:number_groups){
        #inside each block, number_within_each_group
        number_within_each_group <- (n:1)[layer_from_top-1]
        for(i in 1:number_within_each_group){
          #row, each group should start from corresponding prior k'th row,col
          #Note 1st in each group k should add to prior k'th ExitBlock (row,col) nrow of k'th Exitbloxks
          #Else the nrow + ncol of the prior exit block should be added to row of i
          Mat[1,(k-1)*number_within_each_group + i] <- if(i == 1){
            List_Mat_UpperLeftRowCol[[ layer_from_top -1 ]][1,k] + list_DimAllExitBlocks[[ layer_from_top -1]][1,k]
          }else Mat[1,(k-1)*number_within_each_group + i -1] + list_DimAllExitBlocks[[ layer_from_top]][1 , (k-1)*number_within_each_group + i -1]
          #col
          Mat[2,(k-1)*number_within_each_group + i] <- if(i == 1){
            List_Mat_UpperLeftRowCol[[ layer_from_top -1 ]][2,k] + list_DimAllExitBlocks[[ layer_from_top ]][1,(k-1)*number_within_each_group + i]
          } else Mat[2,(k-1)*number_within_each_group + i -1] + list_DimAllExitBlocks[[ layer_from_top]][1,(k-1)*number_within_each_group + i]

        }
      }
    }
    Mat
  }


  #Find upper left corner of placement of ExitBlok (Index, number)

  #List of matrices to hold Row,Col for Exit Matrices. To be filled out by UpdaterFun
  List_Mat_UpperLeftRowCol <- vector('list', length = n)

  #First Index is simply the one large Exit Block, for all variables.
  List_Mat_UpperLeftRowCol[[1]] <- matrix( c(1, prod(unlist(listDim)) + 1),byrow = T,nrow = 2,ncol = 1)

  #For the rest we have that n >= 4 such that we have >= 2 layers.

  for(i in 2:n){
    List_Mat_UpperLeftRowCol[[i]] <- UpdaterFun(list_DimAllExitBlocks = list_DimAllExitBlocks, layer_from_top = i)
  }


  ## Nice to have:
  # List_Mat_UpperLeftRowCol
  # list_DimAllExitBlocks
  # list_IndexAllExitBlocks


  #### Insert ExitMatrices
  a <- BlockGenerator(MatListIndex = seq_len(length(listS)), listS = listS, listAlpha = listAlpha)
  a <- as.matrix(a)

  #Since the smallest block are orderstat of 2 ph, the lowest order Exit block is 3.
  #loop over layers n,...,3
  for(i in seq_len(length(listS)-2) ){
    #looop over every i-length index
    for(k in 1:ncol(List_Mat_UpperLeftRowCol[[i]])){

      block <- ExitBlock(Index = list_IndexAllExitBlocks[[i]][,k])
      dimblock <- list_DimAllExitBlocks[[i]][,k]
      ULrc <- List_Mat_UpperLeftRowCol[[i]][,k]
      a[ULrc[1]:(ULrc[1]+dimblock[1]-1) , ULrc[2]:(ULrc[2]+dimblock[2]-1) ] <- block
    }
  }


  ##### D Exit Matrix ####

  D <- matrix(0, nrow = dim(a)[1] , ncol = factorial(length(listS)) )

  #need row numbers of all 1-exit blocks i.e. exit vectors s.
  for(i in seq_len(ncol(D))){
    exitVec <- ExitBlock(Index = list_IndexAllExitBlocks[[length(listS)]][i])
    rowStart <- List_Mat_UpperLeftRowCol[[length(listS)]][1,i]
    len <- length(exitVec)
    D[rowStart:(rowStart + len -1), i] <- exitVec
  }

  #ALL rows at most has 1 non-zero entry
  if(!all(rowSums(D!=0) %in% c(0,1))) stop("D matrix does not only have 1 entry in each row. or zero")

  ##### pi initial vector ####
  pi <- rep(0, dim(a)[1])
  tmp <- Reduce(f = kronecker, x = listAlpha)
  pi[seq_along(tmp)] <- tmp


  #End of AbsProb fun
  return(list(pi = pi,
              SubMat = a,
              D = D,
              List_Mat_UpperLeftRowCol = List_Mat_UpperLeftRowCol))

}#End of function


#' Dimension of Absorption order matrix
#'
#' @param dimVec vector of dimensions of the phase-types.
#'
#' @return The corresponding dimension of the D matrix for calculating the absorption order probabilities.
#' @export
#'
#' @examples
dimAbsOrderPH <- function(dimVec){

  n <- length(dimVec)

  if(length(dimVec) == 2) {
    sum(dimVec) + prod(dimVec)
  }
  else {
    prod(dimVec) +
      sum(apply(as.matrix(combn(x = dimVec, n-1)), MARGIN = 2, FUN = function(x) dimAbsOrderPH(dimVec =  x) ))
  }
}
