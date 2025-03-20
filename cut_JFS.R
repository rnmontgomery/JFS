# R Function to estimate cut point for JFS

# ARGUMENTS:

# feasible_values: Values for the feasibility criteria that are feasible. If include_recruitment = TRUE, the first value is the include_recruitment rate (monthly) and the second is retention.
# infeasible_values: Values for the feasibility criteria that are infeasible. If include_recruitment = TRUE, the first value is the include_recruitment rate and the second is retention.
# sample_size: Total number of observations in the study.
# preliminary_sample_size: Preliminary sample size if including both monthly include_recruitment rate and retention rate and modeling the trade-offs.
# include_recruitment: TRUE/FALSE, whether include_recruitment rate is included in the analysis.
# recruitment_window: Estimated include_recruitment window for full trial.
# num_simulations: Number of simulated datasets used in the estimation process.
# num_bootstrap: Number of bootstrap replications for uncertainty estimation.
# probability_threshold: Probability threshold for proceeding given the feasibility criteria.
# prior_alpha: Prior prior_alpha parameter for Bayesian estimation.
# prior_beta: Prior prior_beta parameter for Bayesian estimation.

# 

# RETURNS

set.seed(123)
cut_jfs <- function(feasible_values, infeasible_values, sample_size, PpF = 0.80, num_simulations = 1000, prior_alpha = c(1), prior_beta = c(1),include_recruitment = FALSE,
                   recruitment_window = NULL, preliminary_sample_size = NULL,alpha_rec = NULL,beta_rec = NULL) {
    # only able to analyses recruitment if retention is used
    if (max(feasible_values) > 1 & length(feasible_values) == 1) {
        stop("Recruitment needs to be analyses with retention")
    }
    # makes a vector of equal length to feasible_values
    if (length(feasible_values) != length(infeasible_values)) {
        stop("feasible_values and infeasible_values must be the same length.")
    }
    # able to study multiple params at once
    if (length(prior_alpha) == 1) {
        prior_alpha <- rep(prior_alpha, length(feasible_values))
    } else if (length(prior_alpha) != length(feasible_values)) {
        stop("prior_alpha must be either a scalar or the same length as feasible_values.")
    }
    # makes a scalar a vector for prior_alpha
    if (length(prior_beta) == 1) {
        prior_beta <- rep(prior_beta, length(feasible_values))
    } else if (length(prior_beta) != length(feasible_values)) {
        stop("prior_beta must be either a scalar or the same length as feasible_values.")
    }
    
    # only for proportions, only prior_beta-binomial
    if (include_recruitment == FALSE) {
        # error handling that feasible_values
        if (max(feasible_values) > 1 | max(infeasible_values) > 1 | min(feasible_values) < 0 | min(infeasible_values) < 0) {
            stop("Values in feasible_values and infeasible_values must be in [0,1], to calculate the JFS cutpoint for non-binomial data you can edit the function code.")
        }
        # 0 if lower is better and 1 if higher is better
        direction <- ifelse((feasible_values - infeasible_values) > 0, 1, 0)
        # Simulation
        set.seed(321) # will remove later
        all_conditions <- c(feasible_values, infeasible_values)
        condition_probs <- matrix(NA, num_simulations, 2)
        
        # simulate fake data under feasible_values and infeasible_values distribution
        for (i in 1:num_simulations) {
            # Create fake data set
            sim_mat <- matrix(rbinom(sample_size * length(all_conditions), 1, rep(all_conditions, each = sample_size)), 
                              nrow = sample_size, ncol = length(all_conditions))
            # matrix of the postieror proabability of
            marginal_posts <- numeric(length(feasible_values)) # VECTOR VERSION
            for (j in seq_along(feasible_values)){
                posterior_prob <- pbeta(feasible_values[j], prior_alpha[j] + sum(sim_mat[, j]), prior_beta[j] + sample_size - sum(sim_mat[, j]))
                marginal_posts[j] <- ifelse(direction[j] == 1, 1 - posterior_prob, posterior_prob)
            }
            # NEED TO TALK ABOUT WHAT'S HAPPENIGN HERE 
            condition_probs[i, 1] <- prod(marginal_posts)
            marginal_posts_2 <- matrix(NA, length(infeasible_values), 1)
            for (l in (length(feasible_values) + 1):(length(feasible_values) + length(infeasible_values))) {
                set <- (length(feasible_values) + 1):(length(feasible_values) + length(infeasible_values))
                Il <- which(set == l) # getting first col 
                empty_rows <-which(rowSums(is.na(marginal_posts_2)) == ncol(marginal_posts_2))
                if (direction[Il] == 1) {
                    marginal_posts_2[empty_rows[1], 1] <-1 - pbeta(feasible_values[Il],prior_alpha[Il] + sum(sim_mat[, l]),prior_beta[Il] + sample_size - sum(sim_mat[, l]))
                } else if (direction[Il] == 0) {
                    marginal_posts_2[empty_rows[1], 1] <-pbeta(feasible_values[Il],prior_alpha[Il] + sum(sim_mat[, l]),prior_beta[Il] + sample_size - sum(sim_mat[, l]))
                }
            }
            condition_probs[i, 2] <- prod(marginal_posts_2)
        }#num_simulations loop WHAT?
    } else if (include_recruitment == TRUE) {
        if (feasible_values[1] <= 1) {
            warning("The include_recruitment rate is <= 1, confirm this is not another binomial feasibility parameter.")
        }
        if (max(feasible_values[2:length(feasible_values)]) > 1) {
            stop("Values in feasible_values and infeasible_values must be in [0,1], to calculate the JFS cutpoint for non-binomial data you can edit the function code.")
        }
        
        # DO THESE VARIABLES EVEN NEED TO EXIST 
        recF <- feasible_values[1]
        recI <- infeasible_values[1]
        retF <- feasible_values[2]
        retI <- infeasible_values[2]
        
        if (length(feasible_values) > 2) {
            feasible_values <- feasible_values[3:length(feasible_values)]
            infeasible_values <- infeasible_values[3:length(infeasible_values)]
        } else{
            feasible_values <- numeric(0)
            infeasible_values <- numeric(0)
        }
        combinedposts <- matrix(NA, num_simulations, 4)
        print("Estimating joint posterior, this may take a few minutes.") # I THINK WE CAN MAKE THIS FASTER, IF IT'S A VECTORIZED FUNC
        
        # THIS FOR LOOP HAS REPEATED CODE
        # DESIGNING A FUNCITON TO HANDLE THIS MORE CLEANLY 
        for (i in 1:num_simulations) {
            # figure out how long it takes to get n people WRITE A FUNCTION WRITE 1X RUN N 
            # Recruitment feasible_values
            weeklylambda <- recF / 4.33
            recFsim <-rpois(ceiling((sample_size / weeklylambda) * 4), weeklylambda)
            csum <- cumsum(recFsim)
            nl <- which(csum >= sample_size)[1]
            recFsim <- recFsim[1:nl]
            diff <- sum(recFsim) - sample_size
            recFsim[nl] <- recFsim[nl] - diff
            
            posterior_recF <-rgamma(num_simulations,sum(recFsim) + alpha_rec,beta_rec + length(recFsim))
            
            # Recruitment infeasible_values
            weeklylambdaI <- recI / 4.33
            recIsim <-rpois(ceiling((sample_size / weeklylambdaI) * 4), weeklylambdaI)
            csum <- cumsum(recIsim)
            nl <- which(csum >= sample_size)[1]
            recIsim <- recIsim[1:nl]
            diff <- sum(recIsim) - sample_size
            recIsim[nl] <- recIsim[nl] - diff
            
            
            posterior_recI <-rgamma(num_simulations,sum(recIsim) + alpha_rec,beta_rec + length(recIsim))
            
            # Retention feasible_values
            
            retention_simF <- rbinom(sample_size, 1, retF)
            
            posterior_retF <-
                rbeta(num_simulations,sum(retention_simF) + prior_alpha[1],prior_beta[1] + length(retention_simF) - sum(retention_simF))
            
            # Retention infeasible_values
            retention_simI <- rbinom(sample_size, 1, retI)
            
            posterior_retI <-rbeta(num_simulations,sum(retention_simI) + prior_alpha[1],prior_beta[1] + length(retention_simI) - sum(retention_simI))
            
            # Joint posterior feasible_values
            simF <- cbind(posterior_retF, posterior_recF)
            
            retention <- seq(0, 1.0, by = 0.001)
            samp_size <- ceiling(preliminary_sample_size / retention)
            rate <- (samp_size / recruitment_window)
            
            jointposteriorF <- matrix(NA, num_simulations, 1)
            for (l in 1:num_simulations) {
                jointposteriorF[l, 1] <- ifelse(simF[l, 2] * 4.33 >=rate[which(round(retention, 3) == round(simF[l, 1], 3))] , 1, 0)
            }
            
            # Joint posterior infeasible_values
            simI <- cbind(posterior_retI, posterior_recI)
            
            jointposteriorI <- matrix(NA, num_simulations, 1)
            for (p in 1:num_simulations) {
                jointposteriorI[p, 1] <- ifelse(simI[p, 2] * 4.33 >= rate[which(round(retention, 3) == round(simI[p, 1], 3))] , 1, 0)
            }
            
            if (length(feasible_values) == 0) {
                jointposteriorotherF <- NA
                jointposteriorotherI <- NA
                
            } else if (length(feasible_values) > 0) {
                all_conditions <- c(feasible_values, infeasible_values)
                direction <- ifelse((feasible_values - infeasible_values) > 0, 1, 0)
                sim_mat <- matrix(NA, sample_size, length(feasible_values) * 2)
                for (l in 1:length(all_conditions)){
                    sim_mat[, l] <- rbinom(sample_size, 1, all_conditions[l])
                }
                marginalpostsF <- matrix(NA, length(feasible_values), 1)
                for (l in 1:length(feasible_values)) {
                    if (direction[l] == 1) {
                        marginalpostsF[l, 1] <-1 - pbeta(feasible_values[l],prior_alpha[l] + sum(sim_mat[, l]),prior_beta[l] + sample_size - sum(sim_mat[, l]))
                    } else if (direction[l] == 0) {
                        marginalpostsF[l, 1] <-pbeta(feasible_values[l],prior_alpha[l] + sum(sim_mat[, l]),prior_beta[l] + sample_size - sum(sim_mat[, l]))
                    }
                }
                
                jointposteriorotherF <- prod(marginalpostsF)
                
                marginalpostsI <- matrix(NA, length(infeasible_values), 1)
                for (l in (length(feasible_values) + 1):(length(feasible_values) + length(infeasible_values))) {
                    set <- (length(feasible_values) + 1):(length(feasible_values) + length(infeasible_values))
                    Il <- which(set == l)
                    empty_rows <-which(rowSums(is.na(marginalpostsI)) == ncol(marginalpostsI))
                    if (direction[Il] == 1) {
                        marginalpostsI[empty_rows[1], 1] <-
                            1 - pbeta(feasible_values[Il],prior_alpha[Il] + sum(sim_mat[, l]),prior_beta[Il] + sample_size - sum(sim_mat[, l]))
                    } else if (direction[Il] == 0) {
                        marginalpostsI[empty_rows[1], 1] <- pbeta(feasible_values[Il],prior_alpha[Il] + sum(sim_mat[, l]), prior_beta[Il] + sample_size - sum(sim_mat[, l]))
                    }
                }
                
                jointposteriorotherI <- prod(marginalpostsI)
                
            }
            
            # Combine joint posteriors
            
            combinedposts[i, 1] <- mean(jointposteriorF)
            combinedposts[i, 2] <- jointposteriorotherF
            combinedposts[i, 3] <- mean(jointposteriorI)
            combinedposts[i, 4] <- jointposteriorotherI
            #print(i)
        }# nsim loop
        
        xs <- seq(0, 1, 0.001)
        
    }
    if (include_recruitment == TRUE) {
        condition_probs <- matrix(NA, num_simulations, 2)
        if (length(feasible_values) == 0) {
            condition_probs[, 1] <- combinedposts[, 1]
            condition_probs[, 2] <- combinedposts[, 3]
        } else if (length(feasible_values) > 0) {
            condition_probs[, 1] <- combinedposts[, 1] * combinedposts[, 2]
            condition_probs[, 2] <- combinedposts[, 3] * combinedposts[, 4]
        }
    }
    
    breaklist <- c(seq(0, 1, 0.0001))
    FAsimopt_F <- condition_probs[, 1]
    FAsimopt_I <- condition_probs[, 2]
    
    # Check:
    picker <- matrix(NA, 10001, 3)
    picker[, 1] <- seq(0, 1, 0.0001)
    for (l in 1:10001) {
        picker[l, 2] <- length(FAsimopt_F[FAsimopt_F >= picker[l, 1]]) / length(FAsimopt_F)
        picker[l, 3] <- length(FAsimopt_I[FAsimopt_I >= picker[l, 1]]) / length(FAsimopt_I)
    }
    
    cutpoint <- picker[tail(which(picker[, 2] >= PpF), n = 1), 1]
    PF <- length(FAsimopt_F[FAsimopt_F >= cutpoint]) / length(FAsimopt_F)
    PI <- length(FAsimopt_I[FAsimopt_I >= cutpoint]) / length(FAsimopt_I)
    
    if (cutpoint == 0) {
        print("Unable to identify a non-zero cut point.")
    } else {
        ret <- list(cutpoint, PF, PI)
        names(ret) <- c("Cutpoint", "P(P|feasible_values)", "P(P|infeasible_values)")
        return(ret)
    }
}


# Motivating example from paper to match cell [1,7] in Table 3
cut_jfs(
    feasible_values = c(0.8),
    infeasible_values = c(0.7),
    num_simulations = 10000,
    sample_size = 20,
    PpF = 0.80,
    prior_alpha = 1,
    prior_beta = 1
)

# Motivating example from paper matches with cutpoint and probabilities
# ARE THE ARGUMENTS FOR feasible_values AND infeasible_values CORRECT? 
# cut_jfs(
#     feasible_values =  c(5.75, 0.8), 
#     infeasible_values = c(4.44, 0.775),
#     num_simulations = 1000,
#     sample_size = 20,
#     PpF = 0.80,
#     preliminary_sample_size = 165,
#     include_recruitment = TRUE,
#     prior_alpha = c(1),
#     prior_beta = c(1),
#     alpha_rec = 0.01,
#     beta_rec = 0.01,
#     recruitment_window = 36
# )
