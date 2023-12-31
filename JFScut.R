

# R Function to estimate cut point for JFS
# The function as written only works for Binomial data 

# Inputs:

# S: Value (or values) for the endpoints we should proceed for (e.g., 0.80 in paper)
# SN: Value (or values) for endpoints we should not proceed for (e..g, 0.70 in paper)
# N: Sample size
# nsims: Number of simulations (fake data sets)
# nboot: Number of bootstrap replications
# PpS: Probability of proceeding given S
# parallel: Run the function in parallel
# ncores: Number of cores to allocate if running parallel 


# Library foreach and doParallel required to run in parallel
library(foreach)
library(doParallel)


JFScut <- function(S, SN, N, PpS = 0.80, parallel = FALSE, ncores = 2, 
                   nsims = 1000, nboot = 1000){
  
  if(length(S)!=length(SN)){
    stop("S and SN must be the same length.")
  }
  
  if(max(S) > 1 | max(SN) > 1 | min(S) < 0 | min(SN) < 0){
    stop("Values in S and SN must be in [0,1], to calculate the JFS cutpoint for non-binomial data you can edit the function code.")
  }
  
  m <- length(S)
  direction <- rep(ifelse((S-SN) > 0, 1,0),2)
  #direction <- c(direction,direction)
  Sdouble <- c(S,S)
  
  # Simulation
  allcond <- c(S,SN)
  conditionprobs <- matrix(NA,nsims,2)
  if (parallel == FALSE){
  for (i in 1:nsims){
    # Create fake data set
    sim_mat <- matrix(NA,N,m*2)
    for(l in 1:length(allcond))
    {
      sim_mat[,l] <- rbinom(N,1, allcond[l])
    }
    
    # Estimate Probability with bootstrap
    boots <- matrix(NA,nboot,(m*2)+2)
    for(j in 1:nboot){
      resamprows <-sample(1:N,N,replace = TRUE) 
      
      for (k in 1:(m*2)){
        if(direction[k] == 1){
        boots[j,k] <- ifelse(mean(sim_mat[resamprows,k]) >= Sdouble[k],1,0)
        } else if (direction[k] == 0){
          boots[j,k] <- ifelse(mean(sim_mat[resamprows,k]) <= Sdouble[k],1,0)
        }
      }
      boots[j,(m*2)+1] <- ifelse(sum(boots[j,1:length(S)])/length(S) == 1,1,0)
      boots[j,(m*2)+2] <- ifelse(sum(boots[j,(length(S)+1):(2*length(S))])/length(S) == 1,1,0)
      
    }
    conditionprobs[i,1] <- mean(boots[,(m*2)+1])
    conditionprobs[i,2] <- mean(boots[,(m*2)+2])
    
  }#nsims loop
  } else if (parallel == TRUE){
    cl <- makeCluster(ncores, outfile="")
    registerDoParallel(cl )
    res = foreach(i = 1:nsims, .combine = 'cbind') %dopar% {
      
      # Create fake data set
      sim_mat <- matrix(NA,N,m*2)
      for(l in 1:length(allcond))
      {
        sim_mat[,l] <- rbinom(N,1, allcond[l])
      }
      
      # Estimate Probability with bootstrap
      boots <- matrix(NA,nboot,(m*2)+2)
      
      for(j in 1:nboot){
        resamprows <-sample(1:N,N,replace = TRUE) 
        
        
       
        for (k in 1:(m*2)){
          if(direction[k] == 1){
            boots[j,k] <- ifelse(mean(sim_mat[resamprows,k]) >= Sdouble[k],1,0)
          } else if (direction[k] == 0){
            boots[j,k] <- ifelse(mean(sim_mat[resamprows,k]) <= Sdouble[k],1,0)
          }
        }
        boots[j,(m*2)+1] <- ifelse(sum(boots[j,1:length(S)])/length(S) == 1,1,0)
        boots[j,(m*2)+2] <- ifelse(sum(boots[j,(length(S)+1):(2*length(S))])/length(S) == 1,1,0)
        
      }
  
      list(mean(boots[,(m*2)+1]), mean(boots[,(m*2)+2]))
    }
    
    conditionprobs[,1] <- as.numeric(unlist(res[1,]))
    conditionprobs[,2] <- as.numeric(unlist(res[2,]))
    stopCluster(cl)
  }  
  
  picker <- matrix(NA,nboot+1,3)
  picker[,1] <- seq(0,1,1/nboot)
  for (m in 1:(nboot+1)){
    picker[m,2] <- length(conditionprobs[conditionprobs[,1]>=picker[m,1],1])/length(conditionprobs[,1])
    picker[m,3] <- length(conditionprobs[conditionprobs[,2]>=picker[m,1],2])/length(conditionprobs[,2])
  }
  
  cutpoint <- picker[tail(which(picker[,2]>=PpS), n = 1),1]
  PS <- picker[tail(which(picker[,2]>=PpS), n = 1),2]
  PSN <- picker[tail(which(picker[,2]>=PpS), n = 1),3]
  
  if (cutpoint == 0){
    print("Unable to identify a non-zero cut point.")
  } else {
  
  ret <- list(cutpoint, PS, PSN)
  names(ret) <- c("Cutpoint", "P(P|S)", "P(P|SN)")
  return(ret)  
}
}


# Example:
#T his example, takes about 9 minutes to run.
# It results in virtually the same cut point, P(P|S) and P(P|SN) as in the paper for 
# one endpoint
# JFScut(S = c(0.8),
#        SN = c(0.7),
#        nsims = 10000,
#        N = 40, 
#        nboot = 1000, 
#        PpS = 0.84, 
#        parallel = FALSE)

# Cut point to match the result of the bottom right cell of Table 2 in the paper
# N = 119, 3 endpoints, with S = (0.80, 0.80, 0.80), SN = (.70, 0.70, 0.70)
# and P(P|S) = 0.93

# JFScut(S = c(0.8, 0.80, 0.80),
#        SN = c(0.70, 0.70, 0.7),
#        nsims = 10000,
#        N = 119, 
#        nboot = 1000, 
#        PpS = 0.93, 
#        parallel = FALSE)
