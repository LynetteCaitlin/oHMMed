#setwd("~/Documents/Documents - Lynette’s MacBook Pro/Research/Active/MalariaProj/optimized_hmm_scripts")
setwd("/Users/vogl/TeX/ordered_hmm/programsNdata_20221114")
library(ggmcmc)
library(ggplot2)
library(gridExtra)
library(cvms)
library(gtools)
library(mistr)
library(zoo)
library(vcd)

source("MCMC_normal_v7.R")
source("MCMC_poisson_v2.R")


#### general MCMC parameters for all the simulations
iter <- 10000
warmup <- floor(iter / 5) # 20%
thin <- 1
print_params <- FALSE
verbose <- TRUE

###########################
# Example 1
###########################
#  Simulate data: -----------------------------------------------------------------------------------------------------
L3 <- 2^13
true_T3p <- rbind(c(0.99, 0.01, 0),
                 c(0.01, 0.98, 0.01),
                 c(0.0, 0.01, 0.99))

true_betas3 <- 1 / (c(0.2, 1.5, 9))
true_alpha3 <- 1

simdata3full <- hmm_simulate_poisgamma_data(L = L3,
                                            mat_T = true_T3p,
                                            betas = true_betas3,
                                            alpha = true_alpha3)
simdata3 <- simdata3full$data


#Set priors: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----

n2_states_inferred <- 2
prior2_Tp <- generate_random_T(n2_states_inferred)
prior_alpha2 <- true_alpha3 - 0.4
prior_betas2 <- c(2,1)

n3_states_inferred <- 3
prior3_Tp <- generate_random_T(n3_states_inferred)
prior_alpha3 <- true_alpha3 #- 0.4
prior_betas3 <- true_betas3 #+ 0.62

n4_states_inferred <- 4
prior4_Tp <- generate_random_T(n4_states_inferred)
prior_alpha4 <- true_alpha3 - 0.4
prior_betas4 <- c(3,2.5,2,1)

n5_states_inferred <- 5
prior5_Tp <- generate_random_T(n5_states_inferred)
prior_alpha5 <- true_alpha3 - 0.4
prior_betas5 <- c(3.3,2.5,2,1,0.2)

#
#  Run MCMC: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----

res3_p_n2 <- hmm_mcmc_pois(data = simdata3,
                           prior_T = prior2_Tp,
                           prior_betas = prior_betas2,
                           prior_alpha = prior_alpha2,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           init_betas = prior_betas2,
                           init_alpha = prior_alpha2,
                           init_T = prior2_Tp,
                           print_params = print_params,
                           verbose = verbose)

## source("MCMC_poisson_v2.R")
res3_p_n3 <- hmm_mcmc_pois(data = simdata3,
                      prior_T = prior3_Tp,
                      prior_betas = prior_betas3,
                      prior_alpha = prior_alpha3,
                      iter = iter,
                      warmup = warmup,
                      thin = thin,
                      init_betas = prior_betas3+0.4,
                      init_alpha = prior_alpha3+0.5,
                      init_T = prior3_Tp,
                      print_params = print_params,
                      verbose = verbose)

par(mfrow=c(2,2))
plot(res3_p_n3$samples$alpha,type='l')
plot(res3_p_n3$samples$betas[,1],type='l')
plot(res3_p_n3$samples$betas[,2],type='l')
plot(res3_p_n3$samples$betas[,3],type='l')
par(mfrow=c(1,1))

plot(res3_p_n3,simulation = TRUE,true_betas = true_betas3,true_alpha = true_alpha3,true_mat_T = true_T3p,true_states = simdata3full$states)


res3_p_n4 <- hmm_mcmc_pois(data = simdata3,
                           prior_T = prior4_Tp,
                           prior_betas = prior_betas4,
                           prior_alpha = prior_alpha4,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           init_betas = prior_betas4,
                           init_alpha = prior_alpha4,
                           init_T = prior4_Tp,
                           print_params = print_params,
                           verbose = verbose)

res3_p_n5 <- hmm_mcmc_pois(data = simdata3,
                           prior_T = prior5_Tp,
                           prior_betas = prior_betas5,
                           prior_alpha = prior_alpha5,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           init_betas = prior_betas5,
                           init_alpha = prior_alpha5,
                           init_T = prior5_Tp,
                           print_params = print_params,
                           verbose = verbose)





summary(res3_p_n3)
coef(res3_p_n3)
plot(res3_p_n3,simulation = TRUE,true_betas = true_betas3,true_alpha = true_alpha3,true_mat_T = true_T3p,true_states = simdata3full$states)

true_alpha3/true_betas3
as.numeric(coef(res3_p_n3)$alpha)/as.numeric(coef(res3_p_n3)$betas)

x <- conf_mat(N = L3, res = res3_p_n3, plot = TRUE)
plot_confusion_matrix(x$`Confusion Matrix`[[1]],add_sums = TRUE)

plot(res3_p_n4,simulation = TRUE,true_betas = true_betas3,true_alpha = true_alpha3,true_mat_T = true_T3p,true_states = simdata3full$states)


