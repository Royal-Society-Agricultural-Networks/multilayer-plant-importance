
# Functions to calculate stability and to build binary and quantitative jacobian matrices

# Refs:
# Sauve et al. Ecology 
# Neutel et al. 2002 Science, 2007 Nature
# Allesina and Tang 2012 Nature
# Staniczenko et al. 2013 Nature Communications


## stability () ##
# stability() returns the value stab necessary for the system to be stable.
# stab is thus such that the real part of the greatest eigenvalue of J-stab*I = 0

stability <- function(J, s2){
  # J is the community matrix  (i.e. the jacobian matrix, see functions jacobian_binary() or jacobian_quant()) with a null diagonal, dim(J) = S*S
  # s2 is an arbitrary parameter such that the real part of the greatest eigenvalue of J-s2*I is negative
  
  test1 <- (dim(J)[1] == dim(J)[2]) # Is J a square matrix?
  test2 <- FALSE
  if (test1 == TRUE){
    S <- dim(J)[1]
    test2 <- which(diag(J) != vector("numeric", S)) # Does J have a null diagonal?
  }
  
  if ((test1 == TRUE) & (length(test2) == 0)){ # if J is a square matrix with a null diagonal
    S <- dim(J)[1] # S is the number of species in the network
    s1 <- 0
    I <- diag(S)
    E1 <- max(Re(eigen(J-s1*I, only.values = T)$values))
    E2 <- max(Re(eigen(J-s2*I, only.values = T)$values))
    if ((E1>=0) & (E2<0)){ # if s2 is well chosen and the system is not already stable
      while ((s2-s1)>=10^-4){
        stab <-(s1+s2)/2
        E1 <-max(Re(eigen(J-stab*I, only.values = T)$values))
        if (E1>=0){
          s1<-stab
        }
        else {
          s2<-stab
        }
      }
      return(stab)
    }
    
    if (E1<0){
      stop("J corresponds to a stable system.")
    }
    if (E2>=0){
      stop("s2 is not high enough.")
    }
  }
  else { # if J is not a square matrix with a null diagonal
    if (test1 == FALSE){
      stop("J is not a square matrix.")
    }
    if (length(test2) > 0){
      stop("J does not have a null diagonal.")
    }
  }
}

## jacobian_binary() ##
# jacobian_binary() returns the jacobian matrix parametrized as for a binary network

jacobian_binary <- function(m){
  # m is the adjacency matrix of the graph, dim(m) = S*S
  # undirected interactions (or bidirectional) are such that m_ij = m_ji = 1
  # directected interactions (such as antagonistic ones) are such that m_ij = 1 and m_ji = 0 if i feeds on j (j -> i)
  # S being the total number of species
  
  if (dim(m)[1] == dim(m)[2]){ # Is m a square matrix?
    J <- m
    J[which(J < t(J))] <- -t(J)[which(J < t(J))]

    L <- sum(abs(J))/2
    strength <- rnorm(L, sd = 0.1) + 1 # L values for interaction strength magnitude drawn from a normal distribution
    upper_tri_str_index <- which((upper.tri(J) == TRUE) & (J != 0), arr.ind = TRUE) # non-zero elements in the upper triangle matrix
    
    lower_tri_str_index <- cbind(upper_tri_str_index[, 2], upper_tri_str_index[, 1]) # so indices match between upper triangle and lower triangle elements
    J[upper_tri_str_index] <- J[upper_tri_str_index]*strength
    
    strength <- rnorm(L, sd = 0.1) + 1
    J[lower_tri_str_index] <- J[lower_tri_str_index]*strength
    return(J)
  }
  else {
    stop("m is not a square matrix.")
  }  
}

## jacobian_quant() ##
# jacobian_quant() returns the jacobian matrix parametrized with quantitative network data.
# Abundance estimates can be based on Staniczenko et al.'s function GetPreference() that calculates preferences and abundances
# If so, when using jacobian_quant(), please use GetPreference() and cite Staniczenko et al.'s paper as well:
# "The ghost of nestedness in ecological networks"
# by Phillip P.A. Staniczenko, Jason C. Kopp and Stefano Allesina
# Nature Communications 4:1391 doi: 10.1038/ncomms2422 (2013)

jacobian_quant <- function(m, ab){
  # m is the adjacency matrix of the graph, dim(m) = S*S whose non-zero elements are interaction frequencies.
  # undirected interactions (or bidirectional) are such that m_ij = m_ji > 0
  # directected interactions (such as antagonistic ones) are such that m_ij > 0 and m_ji = 0 if i feeds on j (j -> i)
  # S being the total number of species
  # ab is a vector of species abundance. ab might be inferred with Staniczenko et al.'s method (Nat. Commun. 2013)
  
  test1 <- (dim(m)[1] == dim(m)[2]) # Is m a square matrix?
  test2 <- (dim(m)[1] == length(ab)) # Has ab the right length (S)?
  
  if ((test1 == TRUE) & (test2 == TRUE)){ # If m is a square matrix and ab has the right length.
    S <- dim(m)[1]
    J <- m
    J[which(J < t(J))] <- -t(J)[which(J < t(J))]
    
    nonzero_index <- which(J != 0, arr.ind = TRUE)
    c<-0
    for (i in 1:dim(nonzero_index)[1]){
      J[nonzero_index[i, 1], nonzero_index[i, 2]] <- J[nonzero_index[i, 1], nonzero_index[i, 2]]/ab[nonzero_index[i, 2]]
      c<-c+1
    }
    
    return(J)
  }
  if (test1 == FALSE){
    stop("m is not a square matrix.")
  }
  if (test2 == FALSE){
    stop("ab doesn't have the right length.")
  }
}
