#' Function for simulation Markov Chains and Markov Jump Process
#'
#' @param tr
#' @param nJump
#' @param phi0
#' @param RunUntilAbs
#'
#' @return
#' @export
#'
#' @examples
simMC <- function(tr, nJump=10, phi0=NULL, RunUntilAbs = FALSE){

  nStates <-  dim(tr)[1]
  cont    <-  any(diag(tr) < 0)

  if(cont){
    pJump <- matrix(data=0,nrow=nStates,ncol=nStates) #Start filling zeros for no i -> i , for non-abs states
    absStates <- which(diag(tr) == 0)

    for(i in 1:nStates){ if (i %in% absStates) pJump[i,i] <- 1 else pJump[i,-i] <- (-tr[i,-i]/tr[i,i]) }
    #define Jump probabilities from Int-Matrix
  }
  else{pJump<-tr
  absStates <- which(diag(tr) == 1)}

  if(is.null(phi0)) {
    phi0<-c(1,rep(0,nStates-1)) # start in state 1 if not supplied start dist
  }

  states <- rep(0,1) #states during the nJump jumps
  jumps  <- rep(0,1) #Desired number of jumps

  states[1] <- sample(nStates,1,prob=phi0) #starting value, from given initial dist. if NULL start in state 1
  abs <- FALSE
  i <- 1

  #Run the chain
  if( RunUntilAbs == T & sum(absStates) == 0) stop('Cannot run the chain until absorption, since there is no abs. state')

  while(i <= nJump ){
    abs <- (states[i] %in% absStates)  # Has the chain been absorbed?

    if( abs & cont || abs & RunUntilAbs ) { #a cont chain should always break. A discrete chain should only break if we ask it to run Until abs
      break
    }

    if(cont)
    {
      jumps[i+1]<- (-rexp(1)/tr[states[i],states[i]])+ jumps[i]
    }
    #add to previous time, recall   exp(rate = lambda) = exp(1)/lambda
    else{
      jumps[i+1] <- i  #Discrete jumps if Markov Chain
    }

    states[i+1] <- sample(nStates,1,prob=pJump[states[i],]) #Jump prob.
    i <- i + 1
    if(RunUntilAbs == TRUE) nJump <- i
  }

  if(cont){
    return(list(states=states,
                t=jumps,
                abs=abs,
                SubIntMat = tr,
                tau = tail(jumps,1)))
  } else {
    return(list(states=states,
                t=jumps,
                abs=abs,
                Pjmp = pJump,
                tau = tail(jumps,1)))
  }

}

#' SimPH , wrapper for simMC, adding the p+1 state and letting the chain run until absorption
#'
#' @param SubIntMat
#' @param phi0
#'
#' @return
#' @export
#'
#' @examples
simPH <- function(SubIntMat = listS[["S3"]], phi0=NULL){
  #Check if cont or discrete
  cont    <-  any( diag(SubIntMat) < 0)

  if(is.null(phi0)) {
    phi0 <- rep(0,dim(SubIntMat)[1]);phi0[1]<-1
  }

  if(cont){
    Intmat <- cbind(SubIntMat, -rowSums(SubIntMat))
    Intmat <- rbind(Intmat,0)
    IntVec <- c(phi0,0)
  } else { #Discrete
    Intmat <- rbind(SubIntMat,0)
    Intmat <- cbind(Intmat, 1-rowSums(Intmat))
    IntVec <- c(phi0,0)
  }
  if(sum(IntVec) < 1) IntVec[length(IntVec)] <- 1-sum(IntVec)

  simMC(tr = Intmat, phi0 = IntVec, RunUntilAbs = T)
}


#' MPH* simulation
#'
#' @param SubIntMat
#' @param IntVec
#' @param RewardMat
#'
#' @return
#' @export
#'
#' @note
#'
#' Note IntVec and IntMat should be p dimensional, for the underlying chain on $E = 1,2,...,p,p+1$
#'
#' @examples
simMPH <- function(SubIntMat, IntVec = NULL, RewardMat){


  if(is.vector(RewardMat)){
    #Reward Matrix is vector, i.e. Transformation by rewards
    if(dim(SubIntMat)[1] != length(RewardMat)) stop('Check dim of reward and statespace')

  } else{
    if(dim(SubIntMat)[1] != dim(RewardMat)[1]) stop('Check dim of reward and statespace')
  }


  pdim <- dim(SubIntMat)[1]

  ifelse(is.vector(RewardMat),
         ydim <- 1,
         ydim <- dim(RewardMat)[2]
  )

  #For call to simMC
  tr <- cbind(SubIntMat, -rowSums(SubIntMat))
  tr <- rbind(tr,0)

  #If initial vector is null, the SimMC call will assign it as (1,0,...,0)
  if(!is.null(IntVec)){
    IntVec <- c(IntVec,0)
  }

  #SimMC
  a <- simMC(tr = tr, phi0 = IntVec, RunUntilAbs = T)
  y <- rep(NA, ydim)

  for( i in 1:ydim){
    # Time spent in each state 1...p
    OccVec <- sapply(X = 1:pdim, function(x) sum(a$t[which(a$states == x)+1]-a$t[which(a$states == x)]))
    OccVec <- as.numeric(OccVec)

    ifelse(is.vector(RewardMat),
           y[i] <- sum(RewardMat*OccVec),
           y[i] <- sum(RewardMat[,i]*OccVec)
    )
  }
  y
}
#Random (Sub)intensity matrix
#Note IntVec and IntMat should be p dimensional, for the underlying chain on E = 1,2,...,p,p+1

#' Simulates n-outcomes of the mph-class
#'
#' @param n
#' @param SubIntMat
#' @param Phi0
#' @param RewardMat
#'
#' @return
#'
#' @examples
rMPH <- function(n, SubIntMat, Phi0, RewardMat){
  replicate(n = n , expr = simMPH(SubIntMat = SubIntMat, IntVec = Phi0, RewardMat = RewardMat ))
}
