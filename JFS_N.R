# R function to estimate required Sample size

# Requires JFScut function

library(dplyr)

# Inputs:
# n_max: Maximum feasible sample size
# PpF: Probability of proceeding given F
# PpI: Probability of proceeding given I
# F: Values for the feasibility criteria that are feasible
# I: Values for the feasibility criteria that are infeasible
# maxit: Maximum number of iterations (increasing sample size by 1) to check if the original estimate fails
# nsimdata: Number of simulated data sets

JFS_N <- function(n_max, PpF, PpI, F, I, maxit = 2, nsimdata = 20000 ,
                  initalsim = 10000,
                  alpha = 1, beta = 1){

 
  initial_n <- function(n_max, PpF, PpI, F, I){
    
    n <- 10:n_max
    direction <- ifelse((F-I) > 0, 1, 0)
    
    if (direction == 1){
      qs <- qbinom(PpF+0.02,n,F, lower.tail = FALSE)
      options <-  1-pbinom(qs,n,F)
      within <- options[between(options, PpF, 1)]
      pos_vs <- which(options %in% within)
      pos_ns <- n[pos_vs]
      qq <- qs[pos_vs]
      
      qsn <- qbinom(PpI,n,I, lower.tail = FALSE)
      optionsn <-  1-pbinom(qsn,n,I)
      withinsn <- optionsn[between(optionsn, 0, PpI)]
      pos_vsn <- which(optionsn %in% withinsn)
      pos_nsn <- n[pos_vsn]
      qqsn <- qsn[pos_vsn]
      
      # Finding matched quantile, or closest value to match as initial n
      dat <- data.frame(n, qs,qsn,options, optionsn)
      dat$c1 <- ifelse(options >= PpF,0, PpF-options)
      dat$c2 <- ifelse(optionsn <= PpI,0, optionsn-PpI)
      dat$c3 <- abs(qs-qsn)
      dat$err <- dat$c1 + dat$c2 + dat$c3
      
    }else if (direction == 0){
      
      qs <- qbinom(PpF,n,F, lower.tail = TRUE)
      options <-  pbinom(qs,n,F)
      within <- options[between(options, PpF, 1)]
      pos_vs <- which(options %in% within)
      pos_ns <- n[pos_vs]
      qq <- qs[pos_vs]
      
      qsn <- qbinom(PpI,n,I, lower.tail = TRUE)
      optionsn <-  pbinom(qsn,n,I)
      withinsn <- optionsn[between(optionsn, 0, PpI)]
      pos_vsn <- which(optionsn %in% withinsn)
      pos_nsn <- n[pos_vsn]
      qqsn <- qsn[pos_vsn]
      
      # Finding matched quantile, or closest value to match as initial n
      dat <- data.frame(n, qs,qsn,options, optionsn)
      dat$c1 <- ifelse(options >= PpF,0, PpF-options)
      dat$c2 <- ifelse(optionsn <= PpI,0, optionsn-PpI)
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
    m <- length(F)
    diffm <- min(abs(F-I))
    x1 <- initial_n(n_max, PpF^(1/m), PpI^(1/m), F[which(abs(F-I) == diffm)[1]],
                    I[which(abs(F-I) == diffm)[1]])
    n_init <- as.numeric(x1[1])
    n_range <- seq(n_init-15, n_init+15, by = 1)
    n_range <- n_range[n_range>0 & n_range <= n_max]
    print("Getting initial estimates")
    estimated <- matrix(NA,length(n_range),2)
    for (step in 1:length(n_range)){
      
      cut <- JFScut(F ,
                    I ,
                    N = n_range[step],
                    PpF,
                    nsims = initalsim,
                    alpha,
                    beta)
      
      if (length(cut) == 1) { # Checks to make sure we have non-zero cut point
        estimated[step,1] <- NA
        estimated[step,2] <- NA
      } else {
        estimated[step,1] <- n_range[step]
        estimated[step,2] <- as.numeric(cut[3])
      }
    }

    metPpI <- estimated[estimated[,2] <= PpI,]
    n_predicted <- metPpI[1,1]
    
    #Check
    print("Checking Result")
    check <- JFScut(F ,
                    I ,
                    N = n_predicted,
                    PpF = PpF,
                    nsims = nsimdata,
                    alpha,
                    beta)
    
    pass1 <- ifelse(as.numeric(check[2]) >= PpF, TRUE, FALSE)
    pass2 <- ifelse(as.numeric(check[3]) <= PpI, TRUE, FALSE)
    
    passed <- ifelse(pass1 == TRUE & pass2 == TRUE, TRUE, FALSE)
    
    if(passed){
      print("Passed = TRUE")
      # Add some ability to decrease sample size
      go <- ifelse(abs(as.numeric(check[2]) - PpF) < 0.01 & abs(as.numeric(check[3]) - PpI)
                   < 0.01, TRUE, FALSE )
      if (go){
        print("Go = TRUE")
        x <- list(n_predicted, check[2], check[3])
        names(x) <- c("N", "PpF", "PpI")
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
          checkm <- JFScut(F ,
                           I ,
                           N = n_predicted-iter,
                           PpF = PpF,
                           nsims = nsimdata,
                           alpha,
                           beta)
          
          pass1 <- ifelse(as.numeric(checkm[2]) >= PpF, TRUE, FALSE)
          pass2 <- ifelse(as.numeric(checkm[3]) <= PpI, TRUE, FALSE)
          
          passed <- ifelse(pass1 == TRUE & pass2 == TRUE, TRUE, FALSE)
          
          if(passed== TRUE){
            print("Decreased, pass = TRUE")
            go <- ifelse(abs(as.numeric(checkm[2]) - PpF) < 0.01 &
                           abs(as.numeric(checkm[3]) - PpI)< 0.01, TRUE, FALSE)
            if (go){
              print("Decreased, go = TRUE")
              
              n_predicted <- n_predicted-iter
              x <- list(n_predicted, check[2], check[3])
              names(x) <- c("N", "PpF", "PpI")
              break
            } else {
              print("Decreased, go = False")
              
              iter <- iter + 1
              x <- list(n_predicted, check[2], check[3])
              break
            }
            
          } else if(passed == FALSE){
            x <- list(n_predicted, check[2], check[3])
            break
          }
          iter <- iter+1
        } # repeat loop
        
      }
      
    } else {
      # If initial check fails, then iterate through maxit times
      # +1 sample size each time
      
      print("Initial estimate failed, rechecking result")
      iter <- 1
      repeat{
        
        checkm <- JFScut(F ,
                         I ,
                         N = n_predicted+iter,
                         PpF = PpF,
                         nsims = nsimdata,
                         alpha,
                         beta)
        
        pass1 <- ifelse(as.numeric(checkm[2]) >= PpF, TRUE, FALSE)
        pass2 <- ifelse(as.numeric(checkm[3]) <= PpI, TRUE, FALSE)
        
        passed <- ifelse(pass1 == TRUE & pass2 == TRUE, TRUE, FALSE)
        
        if (passed){
          n_predicted <- n_predicted + iter
          x <- list(n_predicted, checkm[2], checkm[3])
          names(x) <- c("N", "PpF", "PpI")
          break
        } 
        iter <- iter + 1
        
        if (iter == maxit){
          print("Failed to find an approximate sample size, try increasing maxit, n_max, or nsimdata and initialsim.")
          x <- list()
          break
        }
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

# Example to match Table 3 value:

# JFS_N(n_max = 200, 
#       PpF = 0.8, 
#       PpI = 0.05,
#       F = c(0.8),
#       I = c(0.7),
#       maxit = 5, 
#       nsimdata = 10000,
#       initalsim = 10000)






