# R Function to estimate cut point for JFS

# ARGUMENTS:

# feasible_values: Values for the feasibility criteria that are feasible. If recruitment = TRUE, the first value is the recruitment rate (monthly) and the second is retention.
# infeasible_values: Values for the feasibility criteria that are infeasible. If recruitment = TRUE, the first value is the recruitment rate and the second is retention.
# sample_size: Total number of observations in the study.
# preliminary_sample_size: Preliminary sample size if including both monthly recruitment rate and retention rate and modeling the trade-offs.
# recruitment: TRUE/FALSE, whether recruitment rate is included in the analysis.
# recruitment_window: Estimated recruitment window for full trial.
# num_simulations: Number of simulated datasets used in the estimation process.
# num_bootstrap: Number of bootstrap replications for uncertainty estimation.
# probability_threshold: Probability threshold for proceeding given the feasibility criteria.
# prior_alpha: Prior prior_alpha parameter for Bayesian estimation.
# prior_beta: Prior prior_beta parameter for Bayesian estimation.

# 

# RETURNS


cut_jfs = function(feasible_values, infeasible_values, sample_size, ppf = 0.80, num_simulations = 1000, prior_alpha = c(1), prior_beta = c(1),
                   recruitment_feasible = NULL, recruitment_infeasible = NULL, retention_feasible = NULL, retention_infeasible = NULL, 
                   recruitment_window = NULL, preliminary_sample_size = NULL, alpha_rec = NULL, beta_rec = NULL) {
    # feasible & infeasible recruitment and retention error handling 
    if (xor(is.null(recruitment_feasible), is.null(recruitment_infeasible))) stop("Both 'recruitment_feasible' and 'recruitment_infeasible' must be provided together.")
    if (xor(is.null(retention_feasible), is.null(retention_infeasible))) stop("Both 'retention_feasible' and 'retention_infeasible' must be provided together.")
    if (!is.null(recruitment_feasible) && !is.numeric(recruitment_feasible)) stop("'recruitment_feasible' must be numeric.")
    if (!is.null(recruitment_infeasible) && !is.numeric(recruitment_infeasible)) stop("'recruitment_infeasible' must be numeric.")
    if (!is.null(retention_feasible) && (!is.numeric(retention_feasible) || retention_feasible <= 0 || retention_feasible >= 1)) stop("'retention_feasible' must be numeric and in (0, 1).")
    if (!is.null(retention_infeasible) && (!is.numeric(retention_infeasible) || retention_infeasible <= 0 || retention_infeasible >= 1)) stop("'retention_infeasible' must be numeric and in (0, 1).")
    # only able to analyses recruitment if retention is used
    if (max(feasible_values) > 1 & length(feasible_values) == 1) stop("Recruitment needs to be analyses with retention")
    # makes a vector of equal length to feasible_values
    if (length(feasible_values) != length(infeasible_values))stop("feasible_values and infeasible_values must be the same length.")
    n = length(feasible_values)
    # prior_alpha & prior_beta can only be scalars or vector length n
    if (length(prior_alpha) == 1) prior_alpha = rep(prior_alpha, n)
    if (length(prior_beta) == 1) prior_beta = rep(prior_beta, n)
    if (length(prior_alpha) != n || length(prior_beta) != n){
        stop("prior_alpha and prior_beta must each be either a scalar or the same length as feasible_values.")
    }
    # error handling that feasible_values
    if (max(feasible_values) > 1 | max(infeasible_values) > 1 | min(feasible_values) < 0 | min(infeasible_values) < 0) {
        stop("Values in feasible_values and infeasible_values must be in [0,1], to calculate the JFS cutpoint for non-binomial data you can edit the function code.")
    }
    
    if (recruitment == FALSE) {
        # ASK DR. MONTGOMERY ON APPENDIDNG RETENTION HERE
        direction = ifelse((feasible_values - infeasible_values) > 0, 1, 0)  # 0 if lower is better and 1 if higher is better
        # Simulation of fake data
        feasible_conditional_probs = numeric(num_simulations)
        infeasible_conditional_probs = numeric(num_simulations)
        for (i in 1:num_simulations) { 
            sim_feasible = rbinom(sample_size, 1, feasible_values)
            sim_infeasible = rbinom(sample_size, 1, infeasible_values)
            marginal_posts = numeric(n)
            marginal_posts_2 = numeric(n)
            
            for (j in seq_along(feasible_values)) {
                post_feasible = pbeta(feasible_values[j], prior_alpha[j] + sum(sim_feasible), prior_beta[j] + sample_size - sum(sim_feasible))
                post_infeasible = pbeta(feasible_values[j], prior_alpha[j] + sum(sim_infeasible), prior_beta[j] + sample_size - sum(sim_infeasible))
                if (direction[j] == 1) {
                    marginal_posts[j] = 1 - post_feasible
                    marginal_posts_2[j] = 1 - post_infeasible
                } else {
                    marginal_posts[j] = post_feasible
                    marginal_posts_2[j] = post_infeasible
                }
            }
            
            feasible_conditional_probs[i] = prod(marginal_posts)
            infeasible_conditional_probs[i] = prod(marginal_posts_2)
        }
    } else { # RECRUITMENT == TRUE
        if (feasible_values[1] <= 1) warning("The recruitment rate is <= 1, confirm this is not another binomial feasibility parameter.")
        if (max(feasible_values[2:n]) > 1) {
            stop("Values in feasible_values and infeasible_values must be in [0,1], 
                 to calculate the JFS cutpoint for non-binomial data you can edit the function code.")
        }

        if (n > 2) { # I'M PRETTY POSITIVE THAT WE SHOULD CUT THIS OUT 
            feasible_values = feasible_values[3:n]
            infeasible_values = infeasible_values[3:n]
        } else {
            feasible_values = numeric(0)
            infeasible_values = numeric(0)
        }
        combinedposts = matrix(NA, num_simulations, 4) # NO NEED FOR A MATRIX 
        print("Estimating joint posterior, this may take a few minutes.") # I THINK WE CAN MAKE THIS FASTER, IF IT'S A VECTORIZED FUNC

        # Recruitment feasible_values
        simulate_recruitment = function(recruitment_rate, sample_size, alpha, beta, num_simulations) {
            weekly_lambda = recruitment_rate / 4.33
            rec_sim = rpois(ceiling((sample_size / weekly_lambda) * 4), weekly_lambda) # is this for a month? 
            nl = which(cumsum(rec_sim) >= sample_size)[1]
            rec_sim = rec_sim[1:nl]
            rec_sim[nl] = rec_sim[nl] - (sum(rec_sim) - sample_size)
            return(rgamma(num_simulations, sum(rec_sim) + alpha, beta + length(rec_sim))) # posterior
        }
        
        for (i in 1:num_simulations) {
            posterior_recF = simulate_recruitment(recruitment_feasible, sample_size, alpha_rec, beta_rec, num_simulations)
            posterior_recI = simulate_recruitment(recruitment_infeasible, sample_size, alpha_rec, beta_rec, num_simulations)

            # Retention feasible_values
            retention_simF = rbinom(sample_size, 1, retF)
            posterior_retF = rbeta(num_simulations,sum(retention_simF) + prior_alpha[1],prior_beta[1] + length(retention_simF) - sum(retention_simF))

            # Retention infeasible_values
            retention_simI = rbinom(sample_size, 1, retI)

                    posterior_retI <-rbeta(num_simulations,sum(retention_simI) + prior_alpha[1],prior_beta[1] + length(retention_simI) - sum(retention_simI))

                    # Joint posterior feasible_values
                    simF = cbind(posterior_retF, posterior_recF)

                    retention = seq(0, 1.0, by = 0.001)
                    samp_size = ceiling(preliminary_sample_size / retention)
                    rate = (samp_size / recruitment_window)

                    jointposteriorF = matrix(NA, num_simulations, 1)
                    for (l in 1:num_simulations) {
                        jointposteriorF[l, 1] = ifelse(simF[l, 2] * 4.33 >=rate[which(round(retention, 3) == round(simF[l, 1], 3))] , 1, 0)
                    }

                    # Joint posterior infeasible_values
                    simI = cbind(posterior_retI, posterior_recI)

                    jointposteriorI = matrix(NA, num_simulations, 1)
                    for (p in 1:num_simulations) {
                        jointposteriorI[p, 1] = ifelse(simI[p, 2] * 4.33 >= rate[which(round(retention, 3) == round(simI[p, 1], 3))] , 1, 0)
                    }

                    if (n == 0) {
                        jointposteriorotherF = NA
                        jointposteriorotherI = NA

                    } else if (n > 0) {
                        all_conditions = c(feasible_values, infeasible_values)
                        direction = ifelse((feasible_values - infeasible_values) > 0, 1, 0)
                        sim_mat = matrix(NA, sample_size, n * 2)
                        for (l in 1:length(all_conditions)){
                            sim_mat[, l] = rbinom(sample_size, 1, all_conditions[l])
                        }
                        marginalpostsF = matrix(NA, n, 1)
                        for (l in 1:n) {
                            if (direction[l] == 1) {
                                marginalpostsF[l, 1] <-1 - pbeta(feasible_values[l],prior_alpha[l] + sum(sim_mat[, l]),prior_beta[l] + sample_size - sum(sim_mat[, l]))
                            } else if (direction[l] == 0) {
                                marginalpostsF[l, 1] <-pbeta(feasible_values[l],prior_alpha[l] + sum(sim_mat[, l]),prior_beta[l] + sample_size - sum(sim_mat[, l]))
                            }
                        }

                        jointposteriorotherF = prod(marginalpostsF)

                        marginalpostsI = matrix(NA, n, 1)
                        for (l in (n + 1):(n + n)) {
                            set = (n + 1):(n + n)
                            Il = which(set == l)
                            empty_rows <-which(rowSums(is.na(marginalpostsI)) == ncol(marginalpostsI))
                            if (direction[Il] == 1) {
                                marginalpostsI[empty_rows[1], 1] <-
                                    1 - pbeta(feasible_values[Il],prior_alpha[Il] + sum(sim_mat[, l]),prior_beta[Il] + sample_size - sum(sim_mat[, l]))
                            } else if (direction[Il] == 0) {
                                marginalpostsI[empty_rows[1], 1] = pbeta(feasible_values[Il],prior_alpha[Il] + sum(sim_mat[, l]), prior_beta[Il] + sample_size - sum(sim_mat[, l]))
                            }
                        }

                        jointposteriorotherI = prod(marginalpostsI)
                    }

                    # Combine joint posteriors

                    combinedposts[i, 1] = mean(jointposteriorF)
                    combinedposts[i, 2] = jointposteriorotherF
                    combinedposts[i, 3] = mean(jointposteriorI)
                    combinedposts[i, 4] = jointposteriorotherI
                    #print(i)
                }# nsim loop

                xs = seq(0, 1, 0.001)

    }
    
    picker_values = seq(0, 1, 0.0001)
    picker_feasible = numeric(length(picker_values))
    picker_infeasible = numeric(length(picker_values))
    
    for (i in seq_along(picker_values)) {
        cutoff = picker_values[i]
        picker_feasible[i] = mean(feasible_conditional_probs >= cutoff)
        picker_infeasible[i] = mean(infeasible_conditional_probs >= cutoff)
    }
    
    cut_index = tail(which(picker_feasible >= ppf), n = 1)
    cut_point = picker_values[cut_index]
    pf = mean(feasible_conditional_probs >= cut_point)
    pi = mean(infeasible_conditional_probs >= cut_point)
    if (cut_point == 0) {
        return("Unable to identify a non-zero cut point.")
    } else {
        return(list("Cut point" = cut_point, "P(P|F)" = pf, "P(P|I)" = pi))
    } 
}


# Motivating example from paper to match cell [1,7] in Table 3
cut_jfs(
    feasible_values = c(0.8),
    infeasible_values = c(0.7),
    num_simulations = 10000,
    sample_size = 20,
    ppf = 0.80,
    prior_alpha = 1,
    prior_beta = 1
)




# UPDATES 
# RETENTION NEEDS TO BE AN OPTIONAL VARIABLE
    # - MAKE IT CLEAR FOR THE USER 
# RECRUITMENT REQUIRED FOR SECOND PART 
# REC I AND REC F NEED TO CHANGE 


# Motivating example from paper matches with cutpoint and probabilities
# ARE THE ARGUMENTS FOR feasible_values AND infeasible_values CORRECT? 
# cut_jfs(
#     feasible_values =  c(5.75, 0.8), 
#     infeasible_values = c(4.44, 0.775),
#     num_simulations = 1000,
#     sample_size = 20,
#     ppf = 0.80,
#     preliminary_sample_size = 165,
#     recruitment = TRUE,
#     prior_alpha = c(1),
#     prior_beta = c(1),
#     alpha_rec = 0.01,
#     beta_rec = 0.01,
#     recruitment_window = 36
# )
