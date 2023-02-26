
## load dependencies:
library(ggmcmc)
library(ggplot2)
library(gridExtra)
library(cvms)
library(gtools)
library(mistr)

## set the source directory 
setwd("~/Documents/GitHub/oHMMed/R")
## load the source file for oHMMed with normal emission densities
source("MCMC_normal.R")



## set general MCMC parameters

iter <- 600               # set number of iterations; note this is redundant since it is equal to the default 
warmup <- floor(iter / 5) # length of burnin is 20% of iter; note this is redundant since it is equal to the default
print_params <- FALSE     # parameters after each iteraion will NOT be printed on the screen
verbose <- TRUE           # progress bar will be shown



###################################
# Appendix Example 1
###################################
# -> add 4 states

#  Simulate 'non- confusing' data: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----
L0 <- 2^13
true_T0 <- rbind(c(0.95, 0.05, 0),
                 c(0.025, 0.95, 0.025),
                 c(0.0, 0.05, 0.95))

true_means0 <- c(-5, 0, 5)
true_sd0  <- 1.5

simdata0full <- hmm_simulate_normal_data(L = L0,
                                         mat_T = true_T0,
                                         means = true_means0,
                                         sigma = true_sd0)
simdata0 <- simdata0full$data

# -> quick view
plot(density(simdata0), main = "")

#  Set priors: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----
n_states_inferred <- 3
prior_T <- generate_random_T(n_states_inferred)
prior_means <- c(-18, -1, 12)
prior_sd <- 3


#  Run MCMC: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----
res0 <- hmm_mcmc_normal(data = simdata0,
                        prior_T = prior_T,
                        prior_means = prior_means,
                        prior_sd = prior_sd,
                        iter = iter,
                        warmup = warmup,
                        thin = thin,
                        print_params = print_params,
                        verbose = verbose)

#  Analyze output: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----
summary(res0)
plot(res0, simulation = TRUE, true_means0, true_sd0, true_T0, simdata0full$states)
# data is stored in res0 object. Plot and summary methods extract the data from the
# "hmm_mcmc_normal" object

#######
## ignore
#######

prior_T
prior_pi=get_pi(prior_T)

joint_T=diag(prior_pi)%*%prior_T
sum(joint_T)

###################################
# Appendix Example 2
###################################


#  Simulate 'confusing' data: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----
L1 <- 2^13
sigma1 <- 0.17
means1 <- c(0.55,1,1.9)
true_T1 <- rbind(c(0.975, 0.025, 0),
                 c(0.016, 0.98, 0.004),
                 c(0.0, 0.015, 0.985))
simdata1full <- hmm_simulate_normal_data(L1,true_T1,means1,sigma1)
simdata1 <- simdata1full$data

# -> quick view
plot(density(simdata1), main="")

#  Set different priors: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----
n2_states_inferred <- 2
prior2_T <- generate_random_T(n2_states_inferred)
prior2_means <- c(0,3)
prior2_sd1 <- 1
prior2_sd2 <- 0.05

n3_states_inferred <- 3
prior3_T <- generate_random_T(n3_states_inferred)
prior3_means <- c(0.1, 0.8, 3)
prior3_sd1 <- 1
prior3_sd2 <- 0.05

#  Run different MCMCs: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----

res_n2_1 <- hmm_mcmc_normal(data = simdata1,
                            prior_T = prior2_T,
                            prior_means = prior2_means,
                            prior_sd = prior2_sd1,
                            iter = iter,
                            warmup = warmup,
                            thin = thin,
                            print_params = print_params,
                            verbose = verbose)


res_n3_1 <- hmm_mcmc_normal(data = simdata1,
                            prior_T = prior3_T,
                            prior_means = prior3_means,
                            prior_sd = prior3_sd1,
                            init_means = means1,
                            init_sd = sigma1,
                            init_T = prior3_T,
                            iter = iter,
                            warmup = warmup,
                            thin = thin,
                            print_params = print_params,
                            verbose = verbose)


res_n2_2 <- hmm_mcmc_normal(data = simdata1,
                            prior_T = prior2_T,
                            prior_means = prior2_means,
                            prior_sd = prior2_sd2,
                            iter = iter,
                            warmup = warmup,
                            thin = thin,
                            print_params = print_params,
                            verbose = verbose)


res_n3_2 <- hmm_mcmc_normal(data = simdata1,
                            prior_T = prior3_T,
                            prior_means = prior3_means,
                            prior_sd = prior3_sd2,
                            iter = iter,
                            warmup = warmup,
                            thin = thin,
                            print_params = print_params,
                            verbose = verbose)

#  Analyze outputs: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----

summary(res_n2_1)
plot(res_n2_1, simulation=TRUE,means1,sigma1,true_T1,simdata1full$states)

summary(res_n2_2)
plot(res_n2_2, simulation=TRUE,means1,sigma1,true_T1,simdata1full$states)

summary(res_n3_1)
plot(res_n3_1, simulation=TRUE,means1,sigma1,true_T1,simdata1full$states)

summary(res_n3_2,simdat1)
plot(res_n3_2, simulation=TRUE,means1,sigma1,true_T1,simdata1full$states)
