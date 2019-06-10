#' Discrete Phase-type representation for first occurrence of k-sub sequence.
#'
#' @param pvec probability vector of unique outcomes
#' @param k_sub_seq sequence to be implemented as DPH. Must start with 1 and increment by one for each new unique outcome.
#'
#' @return Initial vector and sub-transition matrix of the discrete phase-type distribution of the first occurrence.
#' @export
#'
#' @examples
DPHSubSeqGen <- function(pvec = c(0.8, 0.2), k_sub_seq = c(1,1,2,2)){

  #Input checks
  if (any(pvec < 0) ) stop("provided probability vector entries must be non-negative")
  #Check probabilities correspond to k_sub_seq
  if(any(1:length(pvec) != (unique(k_sub_seq))) ) stop("Sub_sequence must be ordered from 1 to d'th unique element")

  #Input cheks over
  #Length of sub-sequence.
  k <- length(k_sub_seq)
  Tmat <- matrix(0, nrow = k, ncol = k)

  #Recall i is when there are i-1 of the needed sub-sequence.
  for(i in k:2){
    Short_cuts <- list()
    list_counter <- 1

    for(m in 1:(i-1)){
      if(all(k_sub_seq[(i-m):(i-1)] == k_sub_seq[1:m])) {
        Short_cuts[[list_counter]] <- k_sub_seq[1:m]
        list_counter <- list_counter + 1
      }
    }
    #Filter out shortcuts that are contained in others.
    #I.e. they must have same append
    Len_Short_cuts <- lapply(X = Short_cuts, FUN = length)
    NextOutcome_Short_cuts <- lapply(X = Len_Short_cuts, FUN = function(x) k_sub_seq[x+1] )

    #Check if the next needed of the short_cuts are equal, only keep the longest.
    Max_indices <- c()
    for(b in unique(unlist(NextOutcome_Short_cuts)) ){
      #Select Short_cuts that have the same next needed.
      Same_index <- which(NextOutcome_Short_cuts == b)
      #Select Next_needed shortcuts with needed value b, and select the longest.
      Max_index <- which(Len_Short_cuts == max(unlist(Len_Short_cuts[Same_index])))
      Max_indices <- c(Max_indices, Max_index )
    }

    #Filtered shurt-cuts
    Filtered_Short_cuts <- Short_cuts[Max_indices]
    Filtered_Len_Short_cuts <- lapply(X = Filtered_Short_cuts, FUN = length)
    Filtered_NextOutcome_Short_cuts <- lapply(X = Filtered_Len_Short_cuts, FUN = function(x) k_sub_seq[x+1] )

    #Insert shortcuts
    for(a in unlist(Filtered_Len_Short_cuts)){

      if(a == (k-1)){
        #Tmat[i,1] <- 1- sum(Tmat[i,]) - pvec[k_sub_seq[a+1] ]
      } else{
        #Go to the next using short cuts, including current length
        Tmat[i, a+2] <- pvec[ k_sub_seq[a+1] ]
      }
    }

    #If no shurtcut goes to first macro state, include this option.
    if(!(k_sub_seq[1] %in% unlist(Filtered_NextOutcome_Short_cuts)) ){
      Tmat[i,2] <- pvec[k_sub_seq[1]]
    }

    #return to 0-macro state, after other jumps have been placed
    if(i == k){
      Tmat[i,1] <- 1 - sum(Tmat[i,]) - pvec[k_sub_seq[k] ]
    } else{
      Tmat[i, 1] <-  1 - sum(Tmat[i,])
    }


  } # end i loop over (k:2)

  #first row simply go to first element of sub sequence.
  Tmat[1,2] <- pvec[k_sub_seq[1]]
  Tmat[1,1] <- 1-pvec[k_sub_seq[1]]

  Tmat
}
