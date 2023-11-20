# R function to estimate required Sample size

# Inputs:

# Requires JFScut function
# n_max: Maximum feasible sample size
# PpS: Probability of proceeding given S
# PpSN: Probability of proceeding given SN
# S: Value (or values) for the endpoints we should proceed for (e.g., 0.80 in paper)
# SN: Value (or values) for endpoints we should not proceed for (e..g, 0.70 in paper)
# maxit: Maximum number of iterations (increasing sample size by 1) to check if the original estimate fails
# bootrep: Number of bootstrap replications for each simulated data set. This should equal the number of bootstrap replications that will be used for the analysis of the data.
# nsimdata: Number of simulated data sets
# paral: Run the function in parallel
# ncor: Number of cores to allocate (if parallel = TRUE) 

# Library dplyr required
library(dplyr)

JFS_N <- function(n_max, PpS, PpSN, S, SN, maxit = 2, bootrep = 1000,
                  paral = FALSE,
                  ncor){
  
  initial_n <- function(n_max, PpS, PpSN, S, SN){
    
    n <- 1:n_max
    direction <- ifelse((S-SN) > 0, 1, 0)
    
    if (direction == 1){
      qs <- qbinom(PpS+0.02,n,S, lower.tail = FALSE)
      options <-  1-pbinom(qs,n,S)
      within <- options[between(options, PpS, 1)]
      pos_vs <- which(options %in% within)
      pos_ns <- n[pos_vs]
      qq <- qs[pos_vs]
      
      qsn <- qbinom(PpSN,n,SN, lower.tail = FALSE)
      optionsn <-  1-pbinom(qsn,n,SN)
      withinsn <- optionsn[between(optionsn, 0, PpSN)]
      pos_vsn <- which(optionsn %in% withinsn)
      pos_nsn <- n[pos_vsn]
      qqsn <- qsn[pos_vsn]
      
      # Finding matched quantile, or closest value to match as initial n
      dat <- data.frame(n, qs,qsn,options, optionsn)
      dat$c1 <- ifelse(options >= PpS,0, PpS-options)
      dat$c2 <- ifelse(optionsn <= PpSN,0, optionsn-PpSN)
      dat$c3 <- abs(qs-qsn)
      dat$err <- dat$c1 + dat$c2 + dat$c3
      
    }else if (direction == 0){ # direction = 1 
      
      qs <- qbinom(PpS,n,S, lower.tail = TRUE)
      options <-  pbinom(qs,n,S)
      within <- options[between(options, PpS, 1)]
      pos_vs <- which(options %in% within)
      pos_ns <- n[pos_vs]
      qq <- qs[pos_vs]
      
      qsn <- qbinom(PpSN,n,SN, lower.tail = TRUE)
      optionsn <-  pbinom(qsn,n,SN)
      withinsn <- optionsn[between(optionsn, 0, PpSN)]
      pos_vsn <- which(optionsn %in% withinsn)
      pos_nsn <- n[pos_vsn]
      qqsn <- qsn[pos_vsn]
      
      # Finding matched quantile, or closest value to match as initial n
      dat <- data.frame(n, qs,qsn,options, optionsn)
      dat$c1 <- ifelse(options >= PpS,0, PpS-options)
      dat$c2 <- ifelse(optionsn <= PpSN,0, optionsn-PpSN)
      dat$c3 <- abs(qs-qsn)
      dat$err <- dat$c1 + dat$c2 + dat$c3
      
    }
    
    matches <- dat[dat$err ==0 ,]
    if(dim(matches)[1] >0){
      chose_n <- as.numeric(matches[which(matches$n == min(matches$n)),1])
    } else {
      chose_n <- dat[dat$err == min(dat$err),1]
    }
    
    return(list(chose_n))
    
  } # End of initial_n
  
  x <- NA
    # First guess by simply taking PpS and PpSN to the 1/mth power
    # and using the function for one sample with smallest difference
    # between S and SN
    m <- length(S)
    diffm <- min(abs(S-SN))
    x1 <- initial_n(n_max, PpS^(1/m), PpSN^(1/m), S[which(abs(S-SN) == diffm)[1]],
                     SN[which(abs(S-SN) == diffm)[1]])
    n_init <- as.numeric(x1[1])
    n_range <- seq(n_init-25, n_init+25, by = 1)
    n_range <- n_range[n_range>0 & n_range <= n_max]
    
    print("Getting initial estimates")
    estimated <- matrix(NA,length(n_range),2)
    for (step in 1:length(n_range)){
      cut <- JFScut(S ,
                    SN ,
                    nsims = 400,
                    N = n_range[step], 
                    nboot = 200, 
                    PpS, 
                    parallel = FALSE)
      if (length(cut) == 1) { # Checks to make sure we have non-zero cut point
        estimated[step,1] <- NA
        estimated[step,2] <- NA
      } else {
      estimated[step,1] <- n_range[step]
      estimated[step,2] <- as.numeric(cut[3])
      }
    }
  
    y <- loess(estimated[,2]~estimated[,1], surface = "direct")
    # plot(estimated[,2]~estimated[,1])
    # j <- order(estimated[,1])
    # lines(n_range[j],y$fitted[j],col="red",lwd=3)
    
    # Need to predict from LOESS curve
    
    predicted_loess <- predict(y, n_range)
    n_predicted <- n_range[which(predicted_loess == predicted_loess[predicted_loess<=PpSN][1])]
    
    #Check
    print("Checking Result")
    check <- JFScut(S ,
                    SN ,
                    nsims = nsimdata,
                    N = n_predicted, 
                    nboot = bootrep, 
                    PpS , 
                    parallel = paral,
                    ncores = ncor)
    
    pass1 <- ifelse(as.numeric(check[2]) >= PpS, TRUE, FALSE)
    pass2 <- ifelse(as.numeric(check[3]) <= PpSN, TRUE, FALSE)
    
    passed <- ifelse(pass1 == TRUE & pass2 == TRUE, TRUE, FALSE)
    if(passed){
      print("Passed = TRUE")
      # Add some ability to decrease sample size
      go <- ifelse(abs(as.numeric(check[2]) - PpS) < 0.015 & abs(as.numeric(check[3]) - PpSN)
                < 0.015, TRUE, FALSE )
      if (go){
        print("Go = TRUE")
        x <- list(n_predicted, check[2], check[3])
        names(x) <- c("N", "PpS", "PpSN")
      } else{
        print("Go = FALSE")
        
        iter <- 1

      repeat{
        
        if (iter == maxit){
           n_predicted <- n_predicted-iter
          x <- list(n_predicted, check[2], check[3])
          break
        }
        print("Checkm, decreasing sample size")
        checkm <- JFScut(S ,
                        SN ,
                        nsims = nsimdata,
                        N = n_predicted-iter, 
                        nboot = bootrep, 
                        PpS , 
                        parallel = paral,
                        ncores = ncor)
        
        pass1 <- ifelse(as.numeric(checkm[2]) >= PpS, TRUE, FALSE)
        pass2 <- ifelse(as.numeric(checkm[3]) <= PpSN, TRUE, FALSE)
        
        passed <- ifelse(pass1 == TRUE & pass2 == TRUE, TRUE, FALSE)
        
        if(passed== TRUE){
          print("Decreased pass = TRUE")
          go <- ifelse(abs(as.numeric(checkm[2]) - PpS) < 0.015 &
                         abs(as.numeric(checkm[3]) - PpSN)< 0.015, TRUE, FALSE)
          if (go){
            print("Decreased go = TRUE")
            
            n_predicted <- n_predicted-iter
            x <- list(n_predicted, check[2], check[3])
            names(x) <- c("N", "PpS", "PpSN")
            break
          } else {
            print("Decreased go = False")
            
            iter <- iter +1
          }
         
        } 
        
      } # repeat loop
      
      }
   
    } else {
      # If initial check fails, this iterate through maxit times
      # +1 sample size each time
      
      print("Initial estimate failed, rechecking result")
      iter <- 1
      repeat{
        
        checkm <- JFScut(S ,
                         SN ,
                         nsims = nsimdata,
                         N = n_predicted + iter, 
                         nboot = bootrep, 
                         PpS, 
                         parallel = paral, 
                         ncores = ncor )
        
        pass1 <- ifelse(as.numeric(checkm[2]) >= PpS, TRUE, FALSE)
        pass2 <- ifelse(as.numeric(checkm[3]) <= PpSN, TRUE, FALSE)
        
        passed <- ifelse(pass1 == TRUE & pass2 == TRUE, TRUE, FALSE)
        
        if (passed){
          n_predicted <- n_predicted + iter
          x <- list(n_predicted, checkm[2], checkm[3])
          names(x) <- c("N", "PpS", "PpSN")
          break
        } 
        if (iter == maxit){
          x <- list()
          break
        }
        iter <- iter + 1
        print(iter)
      } # repeat loop
    } # else loop
  
  if (length(x) >1){
    return(x)
  } else if(length(x)==1){
    print("Failed to find an approximate sample size, try increasing maxit, or n_max")
    return(x)
  }
} # end of JFS_N



# Example:

# From Table 1 we would expect approximately N = 119
 # JFS_N(n_max = 200, 
 #       PpS = 0.8, 
 #       PpSN = 0.05,
 #       S = c(0.8),
 #       SN = c(0.7),
 #       maxit = 5, 
 #       nsimdata = 10000
 #       bootrep = 10000)
