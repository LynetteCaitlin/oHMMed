## load dependencies:
library(ggmcmc)
library(ggplot2)
library(gridExtra)
library(cvms)
library(gtools)
library(mistr)
library(zoo)
library(vcd)
library(cowplot)

## set the source directory 
setwd("~/Documents/GitHub/oHMMed/R")
## load the source file for oHMMed with poisson emission densities
## (note it is dependent on the version with normal emissions)
source("MCMC_normal.R")
source("MCMC_poisson.R")

####################################################################
# Please read the following usage recommendations:
# https://github.com/LynetteCaitlin/oHMMed/blob/main/UsageRecommendations.pdf
# This code is for the example simulations in the above. 
####################################################################

##### Set general MCMC parameters:
iter <- 2000                # set number of iterations; note that the default is 1000
warmup <- floor(iter * 0.4) # length of burnin is 40% of iter; note that the default is 20%
print_params <- FALSE       # parameters after each iteration will NOT be printed on the screen
verbose <- TRUE             # progress bar will be shown, as well as messages

##### Simulate a sequence!:

## set parameters:
L1 <- 2^13  # sequence length
            #transition rate matrix between hidden states:
true_T1 <- rbind(c(0.99, 0.01, 0),
                 c(0.01, 0.98, 0.01),
                 c(0.0, 0.01, 0.99))

true_betas1 <- 1 / (c(0.2, 1.5, 9))  # state specific rate parameters for gamma-poisson emission densities
true_alpha1 <- 1.3                   # shared shape parameters across states for gamma-poisson emission densities

# simulation step:
simdata_full_gamma_poisson <- hmm_simulate_gamma_poisson_data(L = L1,
                                                mat_T = true_T1,
                                                betas = true_betas1,
                                                alpha = true_alpha1)
# then extract the simulated data/observed sequence...:
simdata_gp <- simdata_full_gamma_poisson$data
# ... and have a quick peek at the overall emission density:
hist(simdata_gp, breaks = 50, main = "")


##### Set up inference framework!:

# Set priors and initial values!

# RECOMMENDED PROCEDURE:
#   Initial values and prior values should be set to the same value
#     however, for demonstration purposes we will not set initial values here, so by default, these are then drawn from the priors
#   We here set all the priors according to our recommendations
#   Note that these recommendations regarding priors apply to overall emission densities that resemble overdispersed poisson densities
#         and user discretion is advised for overall emission densities that already resemble mixtures of bell curves

#  all prior betas will likely lie below the 'empirical overall beta', which is:
(sum(simdata_gp) / length(simdata_gp)) / ((mean(simdata_gp)^2) / ((var(simdata_gp) - mean(simdata_gp))))
#  but note that setting them near or lower than the overall mean may be a good idea because of the long tails of the overall distribution:
mean(simdata_gp)

# As per recommendation, we set up a series of inference runs for between 2 and 5 states
#     Remember: We do not know the true number of states!
# The recommendations for the betas are heeded below. 

n2_states_inferred <- 2                           # number of states to be inferred
prior2_T <- generate_random_T(n2_states_inferred) # prior transition matrix, randomly generated
prior_alpha2 <- (mean(simdata_gp)^2) / ((var(simdata_gp) - mean(simdata_gp)) / 2) # recommended prior alpha! changes with number of inferred states!
prior_betas2 <- c(5, 1)

n3_states_inferred <- 3
prior3_T <- generate_random_T(n3_states_inferred)
prior_alpha3 <- (mean(simdata_gp)^2) / ((var(simdata_gp) - mean(simdata_gp)) / 3) # recommended prior alpha! changes with number of inferred states!
prior_betas3 <- c(5,3,1)

n4_states_inferred <- 4
prior4_T <- generate_random_T(n4_states_inferred)
prior_alpha4 <- (mean(simdata_gp)^2) / ((var(simdata_gp) - mean(simdata_gp)) / 4) # recommended prior alpha! changes with number of inferred states!
prior_betas4 <- c(5, 4, 3, 1)

n5_states_inferred <- 5
prior5_T <- generate_random_T(n5_states_inferred)
prior_alpha5 <- (mean(simdata_gp)^2) / ((var(simdata_gp) - mean(simdata_gp)) / 5) # recommended prior alpha! changes with number of inferred states!
prior_betas5 <- c(5, 4, 3, 2, 1)

##### Series of oHMMed inference runs for increasing numbers of states:

res_gp_n2 <- hmm_mcmc_gamma_poisson(data = simdata_gp,
                                  prior_T = prior2_T,
                                  prior_betas = prior_betas2,
                                  prior_alpha = prior_alpha2,
                                  iter = iter,
                                  warmup = warmup,
                                  print_params = print_params,
                                  verbose = verbose)

# Recall: it is recommended to also set: init_betas = prior_betas2; init_alpha = prior_alpha2; init_T = prior2_T


res_gp_n3 <- hmm_mcmc_gamma_poisson(data = simdata_gp,
                                  prior_T = prior3_T,
                                  prior_betas = prior_betas3,
                                  prior_alpha = prior_alpha3,
                                  iter = iter,
                                  warmup = warmup,
                                  print_params = print_params,
                                  verbose = verbose)

res_gp_n4 <- hmm_mcmc_gamma_poisson(data = simdata_gp,
                                  prior_T = prior4_T,
                                  prior_betas = prior_betas4,
                                  prior_alpha = prior_alpha4,
                                  iter = iter,
                                  warmup = warmup,
                                  print_params = print_params,
                                  verbose = verbose)

res_gp_n5 <- hmm_mcmc_gamma_poisson(data = simdata_gp,
                                  prior_T = prior5_T,
                                  prior_betas = prior_betas5,
                                  prior_alpha = prior_alpha5,
                                  iter = iter,
                                  warmup = warmup,
                                  print_params = print_params,
                                  verbose = verbose)


#### Check the results:

# first compare the log.likelihoods for different numbers of states
#   and find the plateau
vll1_gp <- c(mean(res_gp_n2$estimates$log_likelihood),
             mean(res_gp_n3$estimates$log_likelihood),
             mean(res_gp_n4$estimates$log_likelihood),
             mean(res_gp_n5$estimates$log_likelihood))

vll2_gp <- c(median(res_gp_n2$estimates$log_likelihood),
             median(res_gp_n3$estimates$log_likelihood),
             median(res_gp_n4$estimates$log_likelihood),
             median(res_gp_n5$estimates$log_likelihood))

ind_state <- c(2, 3, 4, 5)
mat_liks1 <- data.frame(cbind(ind_state, vll1_gp, vll2_gp))

p1 <- ggplot(mat_liks1) + 
  geom_point(aes(ind_state, vll1_gp), size = 3) +
  geom_point(aes(ind_state, vll2_gp), size = 1.5, col = "grey") + 
  geom_hline(yintercept = -15670.33, linetype = "dashed", colour = "grey") + 
  xlab("Number of States") + 
  ylab("Mean (Median) Log-Likelihood")
p1



# looks like the optimal number of states is 3 (where the plateau starts), so examine the results:
#   the summary contains all the estimates, and the approximate Kullback-Leibler divergence
summary(res_gp_n3)
#   graphical diagnostics and confusion matrix
plot(res_gp_n3, simulation = TRUE, true_alpha1, 
     true_betas1, true_T1, simdata_full_gamma_poisson$states)

