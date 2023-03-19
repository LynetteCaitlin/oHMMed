# Saving example data
library(oHMMed)

# 1) Normal model: -------------------------------------------------------------
iter <- 1500          
warmup <- floor(iter * 0.4)
print_params <- FALSE     
verbose <- TRUE
L1 <- 2^13                      
true_sigma1 <- 0.195            
true_means1 <- c(0.55, 1, 1.9)  
true_T1 <- rbind(c(0.7, 0.3, 0),    
                 c(0.3, 0.6, 0.1),
                 c(0.0, 0.35, 0.65))

simdata1full <- hmm_simulate_normal_data(L1, true_T1, true_means1, true_sigma1)

simdata1 <- simdata1full$data 
plot(density(simdata1), main = "")
pr_means <- c(0.6, 0.95, 1.95)

n3_states_inferred <- 3
rand_T <- generate_random_T(n3_states_inferred)


example_hmm_mcmc_normal <- hmm_mcmc_normal(data = simdata1,
                                           prior_T = rand_T,
                                           prior_means = pr_means,
                                           prior_sd = sd(simdata1) / 3,
                                           init_T = rand_T,
                                           init_means = pr_means,
                                           init_sd = sd(simdata1) / 3,
                                           iter = iter,
                                           warmup = warmup,
                                           print_params = print_params,
                                           verbose = verbose)


# 2) Gamma-Poisson model: ------------------------------------------------------

iter <- 1500              
warmup <- floor(iter * 0.4)
print_params <- FALSE     
verbose <- TRUE           
L1 <- 2^13  # sequence length
true_T1 <- rbind(c(0.99, 0.01, 0),
                 c(0.01, 0.98, 0.01),
                 c(0.0, 0.01, 0.99))

true_betas1 <- 1 / (c(0.2, 1.5, 9)) 
true_alpha1 <- 1.3                  

# simulation step:
simdata1full <- hmm_simulate_gamma_poisson_data(L = L1,
                                                mat_T = true_T1,
                                                betas = true_betas1,
                                                alpha = true_alpha1)

simdata1 <- simdata1full$data
hist(simdata1, breaks = 50, main = "")

n3_states_inferred <- 3
prior3_T <- generate_random_T(n3_states_inferred)
prior_alpha3 <- (mean(simdata1)^2) / ((var(simdata1) - mean(simdata1)) / 3) # recommended prior alpha! changes with number of inferred states!
prior_betas3 <- c(5,3,1)

example_hmm_mcmc_gamma_poisson <- hmm_mcmc_gamma_poisson(data = simdata1,
                                                         prior_T = prior3_T,
                                                         prior_betas = prior_betas3,
                                                         prior_alpha = prior_alpha3,
                                                         iter = iter,
                                                         warmup = warmup,
                                                         print_params = print_params,
                                                         verbose = verbose)


# 3) Save examples: ------------------------------------------------------------
# save(list = c("example_hmm_mcmc_normal", "example_hmm_mcmc_gamma_poisson"),
#      file = "data/simulated_models_examples.RData")

