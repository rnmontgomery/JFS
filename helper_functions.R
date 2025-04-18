# HELPER FUNCTIONS: 

# Recruitment feasible_values
simulate_recruitment = function(recruitment_rate, sample_size, alpha, beta, num_simulations) {
    weekly_lambda = recruitment_rate / 4.33
    rec_sim = rpois(ceiling((sample_size / weekly_lambda) * 4), weekly_lambda) # is this for a month? 
    nl = which(cumsum(rec_sim) >= sample_size)[1]
    rec_sim = rec_sim[1:nl]
    rec_sim[nl] = rec_sim[nl] - (sum(rec_sim) - sample_size)
    return(rgamma(num_simulations, sum(rec_sim) + alpha, beta + length(rec_sim))) # posterior
}


############################################################################################################
# NOTE: THE get_conditional_probs FUNCTION is expecting vectors of priors that the error handling will fix 
############################################################################################################

# OG version after reducing from using matrices 
get_conditional_probs = function(n, feasible_values, infeasible_values, sample_size, prior_alpha, prior_beta, num_simulations) {
    direction = ifelse((feasible_values - infeasible_values) > 0, 1, 0)
    feasible_conditional_probs = numeric(num_simulations)
    infeasible_conditional_probs = numeric(num_simulations)
    
    for (i in 1:num_simulations) {
        sim_feasible = rbinom(sample_size, 1, feasible_values)
        sim_infeasible = rbinom(sample_size, 1, infeasible_values)
        
        post_feasible = pbeta(feasible_values, prior_alpha + sum(sim_feasible), prior_beta + sample_size - sum(sim_feasible))
        post_infeasible = pbeta(feasible_values, prior_alpha + sum(sim_infeasible), prior_beta + sample_size - sum(sim_infeasible))
        
        marginal_posts_f = ifelse(direction == 1, 1 - post_feasible, post_feasible)
        marginal_posts_i = ifelse(direction == 1, 1 - post_infeasible, post_infeasible)
        
        feasible_conditional_probs[i] = prod(marginal_posts_f)
        infeasible_conditional_probs[i] = prod(marginal_posts_i)
    }
    
    return(list(
        feasible_conditional_probs = feasible_conditional_probs,
        infeasible_conditional_probs = infeasible_conditional_probs
    ))
}

# new vectorized version 
get_conditional_probs_vecorized = function(feasible_values, infeasible_values, n, sample_size, prior_alpha, prior_beta,num_simulations) {
    direction = ifelse((feasible_values - infeasible_values) > 0, 1, 0)
    
    simulate_once = function() {
        sim_feasible = rbinom(sample_size, 1, feasible_values)
        sim_infeasible = rbinom(sample_size, 1, infeasible_values)
        
        post_feasible = pbeta(
            feasible_values,
            prior_alpha + sum(sim_feasible),
            prior_beta + sample_size - sum(sim_feasible)
        )
        
        post_infeasible = pbeta(
            feasible_values,
            prior_alpha + sum(sim_infeasible),
            prior_beta + sample_size - sum(sim_infeasible)
        )
        
        marginal_posts_f = ifelse(direction == 1, 1 - post_feasible, post_feasible)
        marginal_posts_i = ifelse(direction == 1, 1 - post_infeasible, post_infeasible)
        
        return(c(
            feasible = prod(marginal_posts_f),
            infeasible = prod(marginal_posts_i)
        ))
    }
    
    results = sapply(1:num_simulations, function(x) simulate_once())
    
    return(list(
        feasible_conditional_probs = results["feasible", ],
        infeasible_conditional_probs = results["infeasible", ]
    ))
}




