# R Function to estimate cut point for JFS

# Inputs:

# F: Values for the feasibility criteria that are feasible, if recruitment = TRUE the first value is the recruitment rate (monthly) and the second is retention
# I: Values for the feasibility criteria that are infeasible, if recruitment = TRUE the first value is the recruitment rate and the second is retention
# N: Sample size
# nprelim: Preliminary sample size if including both monthly recruitment rate and retention rate and modeling the trade-offs
# recruitment: TRUE/FALSE, is recruitment rate included 
# nprelim: Preliminary sample size estimate
# recwindow: Estimated recruitment window for full trial
# nsims: Number of simulations (fake data sets)
# nboot: Number of bootstrap replications
# Ppf: Probability of proceeding given F
# alpha: Prior alpha parameter
# beta: Prior beta parameter



library(dplyr)
library(Hmisc)

JFScut <- function(F, I, N, PpF = 0.80, 
                   nsims = 1000, alpha = c(1), beta = c(1),
                   recruitment = FALSE, recwindow = NULL,
                   nprelim = NULL, alpha_rec = NULL, beta_rec = NULL){
  
  if(max(F) > 1 & length(F) == 1){
    stop("Recruitment needs to be analyses with retention")
  }
  
  if(length(F)!=length(I)){
    stop("F and I must be the same length.")
  }
  
  if (length(F)>length(alpha)){
    alpha <- rep(alpha, length(F))
  }
  
  if (length(F)>length(beta)){
    beta <- rep(beta, length(F))
  }
  
  if (recruitment == FALSE){
    
    if(max(F) > 1 | max(I) > 1 | min(F) < 0 | min(I) < 0){
      stop("Values in F and I must be in [0,1], to calculate the JFS cutpoint for non-binomial data you can edit the function code.")
    }
    
    if (F[1]>1){
      stop("The first value in F is > 1, should recruitment be set to TRUE?")
    }
    
    m <- length(F)
    direction <- ifelse((F-I) > 0, 1, 0)
    
    # Simulation
    allcond <- c(F,I)
    conditionprobs <- matrix(NA,nsims,2)
    
    for (i in 1:nsims){
      # Create fake data set
      sim_mat <- matrix(NA,N,m*2)
      for(l in 1:length(allcond))
      {
        sim_mat[,l] <- rbinom(N,1, allcond[l])
      }
      
      marginalposts <- matrix(NA,length(F),1)
      for (l in 1:length(F)){
        
        if (direction[l] == 1){
          marginalposts[l,1] <- 1-pbeta(F[l],alpha[l] + sum(sim_mat[,l]), beta[l] + N - sum(sim_mat[,l]))
        } else if (direction[l] == 0){
          marginalposts[l,1] <- pbeta(F[l],alpha[l] + sum(sim_mat[,l]), beta[l] + N - sum(sim_mat[,l]))
        }
      }
      
      conditionprobs[i,1] <-prod(marginalposts)
      marginalposts2 <- matrix(NA,length(I),1)
      for (l in (length(F)+1):(length(F)+ length(I))){
        set <- (length(F)+1):(length(F)+ length(I))
        Il <- which(set == l)
        empty_rows <- which(rowSums(is.na(marginalposts2)) == ncol(marginalposts2))
        if (direction[Il] == 1){
          marginalposts2[empty_rows[1],1] <- 1-pbeta(F[Il],alpha[Il] + sum(sim_mat[,l]), beta[Il] + N - sum(sim_mat[,l]))
        } else if(direction[Il] == 0){
          marginalposts2[empty_rows[1],1] <- pbeta(F[Il],alpha[Il] + sum(sim_mat[,l]), beta[Il] + N - sum(sim_mat[,l]))
        }
      }
      
      conditionprobs[i,2] <-prod(marginalposts2)
      
      # prod(marginalposts)
      # prod(marginalposts2)
      # colMeans(sim_mat)
      
    }#nsims loop
  } else if (recruitment == TRUE){
    
    if (F[1]<=1){
      warning("The recruitment rate is <= 1, confirm this is not another binomial feasibility parameter.")
    }
    
    if (max(F[2:length(F)])>1){
      stop("Values in F and I must be in [0,1], to calculate the JFS cutpoint for non-binomial data you can edit the function code.")
    }
    
    recF <- F[1]
    recI <- I[1]
    retF <- F[2]
    retI <- I[2]
    
    if (length(F)>2){
    F <- F[3:length(F)]
    I <- I[3:length(I)]
    } else{
      F <- numeric(0)
      I <- numeric(0)
    }
    combinedposts <- matrix(NA,nsims,4)
    print("Estimating joint posterior, this may take a few minutes.")
    
    for (i in 1:nsims){
      
      # Recruitment F
      weeklylambda <- recF/4.33
      recFsim <- rpois(ceiling((N/weeklylambda)*4),weeklylambda )
      csum <- cumsum(recFsim)
      nl <- which(csum>=N)[1]
      recFsim <- recFsim[1:nl]
      diff <- sum(recFsim)-N
      recFsim[nl] <- recFsim[nl]-diff
      
      posterior_recF <- rgamma(nsims, sum(recFsim) + alpha_rec, 
                               beta_rec + length(recFsim))
      
      # Recruitment I
      weeklylambdaI <- recI/4.33
      recIsim <- rpois(ceiling((N/weeklylambdaI)*4),weeklylambdaI )
      csum <- cumsum(recIsim)
      nl <- which(csum>=N)[1]
      recIsim <- recIsim[1:nl]
      diff <- sum(recIsim)-N
      recIsim[nl] <- recIsim[nl]-diff
      
      
      posterior_recI <- rgamma(nsims, sum(recIsim) + alpha_rec, 
                               beta_rec + length(recIsim))
      
      # Retention F
      
      retention_simF <- rbinom(N,1,retF)
      
      posterior_retF <- rbeta(nsims, sum(retention_simF) + alpha[1], 
                              beta[1] + length(retention_simF) - sum(retention_simF))
      
      
      # Retention I
      retention_simI <- rbinom(N,1,retI)
      
      posterior_retI <- rbeta(nsims, sum(retention_simI) + alpha[1], 
                              beta[1] + length(retention_simI) - sum(retention_simI))
      
      
      # Joint posterior F
      simF <- cbind(posterior_retF, posterior_recF)
      
      retention <- seq(0, 1.0, by = 0.001)
      samp_size <- ceiling(nprelim/retention)
      rate <- (samp_size/recwindow)
      
      jointposteriorF <- matrix(NA,nsims,1)
      for (l in 1:nsims){
        jointposteriorF[l,1] <- ifelse(simF[l,2]*4.33 >= 
                                         rate[which(round(retention,3) == round(simF[l,1],3))] ,1,0)
      }
      
      
      # Joint posterior I
      simI <- cbind(posterior_retI, posterior_recI)
      
      
      jointposteriorI <- matrix(NA,nsims,1)
      for (p in 1:nsims){
        jointposteriorI[p,1] <- ifelse(simI[p,2]*4.33 >= 
                                         rate[which(round(retention,3) == round(simI[p,1],3))] ,1,0)
      }
      
      if (length(F) == 0){
        
        jointposteriorotherF <- NA
        jointposteriorotherI <- NA
        
      } else if (length(F) > 0){
        
        
        allcond <- c(F,I)
        direction <- ifelse((F-I) > 0, 1, 0)
        
        
        sim_mat <- matrix(NA,N,length(F)*2)
        for(l in 1:length(allcond))
        {
          sim_mat[,l] <- rbinom(N,1, allcond[l])
        }
        
        marginalpostsF <- matrix(NA,length(F),1)
        for (l in 1:length(F)){
          
          if (direction[l] == 1){
            marginalpostsF[l,1] <- 1-pbeta(F[l],alpha[l] + sum(sim_mat[,l]), beta[l] + N - sum(sim_mat[,l]))
          } else if (direction[l] == 0){
            marginalpostsF[l,1] <- pbeta(F[l],alpha[l] + sum(sim_mat[,l]), beta[l] + N - sum(sim_mat[,l]))
          }
        }
        
        jointposteriorotherF <-prod(marginalpostsF)
        
        marginalpostsI <- matrix(NA,length(I),1)
        for (l in (length(F)+1):(length(F)+ length(I))){
          set <- (length(F)+1):(length(F)+ length(I))
          Il <- which(set == l)
          empty_rows <- which(rowSums(is.na(marginalpostsI)) == ncol(marginalpostsI))
          if (direction[Il] == 1){
            marginalpostsI[empty_rows[1],1] <- 1-pbeta(F[Il],alpha[Il] + sum(sim_mat[,l]), beta[Il] + N - sum(sim_mat[,l]))
          } else if(direction[Il] == 0){
            marginalpostsI[empty_rows[1],1] <- pbeta(F[Il],alpha[Il] + sum(sim_mat[,l]), beta[Il] + N - sum(sim_mat[,l]))
          }
        }
        
        jointposteriorotherI <-prod(marginalpostsI)
        
      }
      
      # Combine joint posteriors
      
      combinedposts[i,1] <- mean(jointposteriorF)
      combinedposts[i,2] <- jointposteriorotherF
      combinedposts[i,3] <- mean(jointposteriorI)
      combinedposts[i,4] <- jointposteriorotherI
      #print(i)
    }# nsim loop
    
    xs <- seq(0,1,0.001)
    
    
  }
  if (recruitment == TRUE){
    conditionprobs <- matrix(NA,nsims,2)
    if (length(F) == 0){
      conditionprobs[,1] <- combinedposts[,1]
      conditionprobs[,2] <- combinedposts[,3]
    } else if(length(F) > 0){
      conditionprobs[,1] <- combinedposts[,1]*combinedposts[,2]
      conditionprobs[,2] <- combinedposts[,3]*combinedposts[,4]
    }
  }
  
  breaklist <- c(seq(0,1, 0.0001))
  FAsimopt_F <- conditionprobs[,1]
  FAsimopt_I <- conditionprobs[,2]
  
  
  # Check:
  picker <- matrix(NA,10001,3)
  picker[,1] <- seq(0,1,0.0001)
  for (l in 1:10001){
    
    picker[l,2] <- length(FAsimopt_F[FAsimopt_F>=picker[l,1]])/length(FAsimopt_F)
    picker[l,3] <- length(FAsimopt_I[FAsimopt_I>=picker[l,1]])/length(FAsimopt_I)
  }
  
  cutpoint <- picker[tail(which(picker[,2]>=PpF), n = 1),1]
  PF <- length(FAsimopt_F[FAsimopt_F>=cutpoint])/length(FAsimopt_F) 
  PI <- length(FAsimopt_I[FAsimopt_I>=cutpoint])/length(FAsimopt_I)
  
  if (cutpoint == 0){
    print("Unable to identify a non-zero cut point.")
  } else {
    ret <- list(cutpoint, PF, PI)
    names(ret) <- c("Cutpoint", "P(P|F)", "P(P|I)")
    return(ret)  
  }
}


# Motivating example from paper to match cell [1,7] in Table 3
# JFScut(F = c(0.8),
#        I = c(0.7),
#        nsims = 10000,
#        N = 20, 
#        PpF = 0.80,
#        alpha = 1,
#        beta = 1)
# 
# # Motivating example from paper matches with cutpoint and probabilities
# JFScut(F =  c(5.75, 0.8),
#        I = c(4.44, 0.775),
#        nsims = 1000,
#        N = 20, 
#        PpF = 0.80,
#        nprelim = 165,
#        recruitment = TRUE,
#        alpha = c(1),
#        beta = c(1),
#        alpha_rec = 0.01,
#        beta_rec = 0.01,
#        recwindow = 36)
