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

####################################################################
# Please read the following usage recommendations:
# https://github.com/LynetteCaitlin/oHMMed/UsageRecommendations.pdf
# This code is for the example simulations in the above. 
####################################################################

##### Set general MCMC parameters:
iter <- 600               # set number of iterations; note this is redundant since it is equal to the default 
warmup <- floor(iter / 5) # length of burnin is 20% of iter; note this is redundant since it is equal to the default
print_params <- FALSE     # parameters after each iteration will NOT be printed on the screen
verbose <- TRUE           # progress bar will be shown


##### Simulate a sequence!:

## set parameters:
L1=2^13                   # sequence length
true_sigma1 <- 0.195           # standard deviation for emission densities of all states
true_means1 <- c(0.55,1,1.9)   # state-specific means for emission densities
      #transition rate matrix between hidden states:
true_T1 <- rbind(c(0.7, 0.3, 0),    
                 c(0.3, 0.6, 0.1),
                 c(0.0, 0.35, 0.65))

# simulation step:
simdata1full <- hmm_simulate_normal_data(L1,true_T1,true_means1,true_sigma1)
# then extract the simulated data/observed sequence...:
simdata1 <- simdata1full$data 
# ... and have a quick peek at the overall emission density:
plot(density(simdata1), main="")


##### Set up inference framework!:

# Set priors and initial values!

# RECOMMENDED PROCEDURE:
#   initial values and prior values should be set to the same value

# run:
c(min(simdata1),max(simdata1))
# then, set initial and prior means to either:
#     (1) values near modes or saddle points seen in the previous density plot
#  or (2) equally spaced values within the just calculated range of the data
# for example (using (1)) and assuming 3 states:
pr_means <- c(0.6,0.95,1.95)

# run:
sd(simdata1)
# then, set the initial and prior standard deviation to the above value divided by the number of states

# the prior/initial transition rate matrix can be randomly generated
# for eg, by assuming 3 states:
rand_T <- generate_random_T(n3_states_inferred)

# overall, the recommended procedure would lead us to run the following inference procedure for 3 states:
res_opt_n3 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = rand_T,
                           prior_means = pr_means,
                           prior_sd = sd(simdata1)/3,
                           init_T = rand_T,
                           init_means = pr_means,
                           init_sd = sd(simdata1)/3,
                           iter = iter,
                           warmup = warmup,
                           print_params = print_params,
                           verbose = verbose)

# note that this should be repeated with different numbers of states, for eg. 2,3,4,5 (see below)
# here, we know that we have simulated 3 hidden states, but in general we would not know the optimal number

# check the results:
#   the summary contains all the estimates, and the approximate kullback-leibler divergence
summary(res_opt_n3)
#   graphical diagnostics and confusion matrix
plot(res_opt_n3, simulation=TRUE,true_means1,true_sigma1,true_T1,simdata1full$states)

# Note:
# if this were not a simulation, but a real example, the graphical diagnostics are: 
#     plot(res_opt_n3)
# and use the following confusion matrix to assess stability of parameter ranges:
#    x <- conf_mat(N=L1,res=opt_n3)
#    plot_confusion_matrix(x$`Confusion Matrix`[[1]],add_sums = TRUE)



# DEMONSTRATION:
# For demonstration purposes, we use different priors:
# As per recommendation, we set up a series of inference runs for between 2 and 5 states
#     Remember: We do not know the true number of states!
# But we will not set the initial values; so by default, they will be sampled from the prior distribution
# Further, we will set a higher prior standard deviation than recommended and more poorly placed means
# So, the prior parameters are set as the following per number of states:

n2_states_inferred <- 2                          # number of states to be inferred
prior2_T <- generate_random_T(n2_states_inferred)# prior transition matrix, randomly generated
prior2_means <- c(0,3)                           # prior means
prior2_sd <- 2.5/2                               # prior standard deviation

n3_states_inferred <- 3
prior3_T <- generate_random_T(n3_states_inferred)
prior3_means <- c(0.1, 0.8, 3)
prior3_sd <- 2.5/3

n4_states_inferred <- 4
prior4_T <- generate_random_T(n4_states_inferred)
prior4_means <- c(0.1,0.7, 0.8, 3)
prior4_sd <-2.5/4

n5_states_inferred <- 5
prior5_T <- generate_random_T(n5_states_inferred)
prior5_means <- c(0.1,0.7, 0.8, 1.5, 3)
prior5_sd <- 2.5/5


##### Series of oHMMed inference runs for increasing numbers of states:

res1_n2 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior2_T,
                           prior_means = prior2_means,
                           prior_sd = prior2_sd,
                           iter = iter,
                           warmup = warmup,
                           print_params = print_params,
                           verbose = verbose)

res1_n3 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior3_T,
                           prior_means = prior3_means,
                           prior_sd = prior3_sd,
                           warmup = warmup,
                           print_params = print_params,
                           verbose = verbose)

res1_n4 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior4_T,
                           prior_means = prior4_means,
                           prior_sd = prior4_sd,
                           warmup = warmup,
                           print_params = print_params,
                           verbose = verbose)

res1_n5 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior5_T,
                           prior_means = prior5_means,
                           prior_sd = prior5_sd,
                           iter = iter,
                           warmup = warmup,
                           print_params = print_params,
                           verbose = verbose)


##### Check the results:

# first compare the log.likelihoods for different numbers of states
#   and find the plateau
vll1=c(mean(res1_n2$estimates$log_likelihood),
       mean(res1_n3$estimates$log_likelihood),
       mean(res1_n4$estimates$log_likelihood),
       mean(res1_n5$estimates$log_likelihood))

vll2=c(median(res1_n2$estimates$log_likelihood),
       median(res1_n3$estimates$log_likelihood),
       median(res1_n4$estimates$log_likelihood),
       median(res1_n5$estimates$log_likelihood))

ind <- c(2,3,4,5)
mat_liks1 <- data.frame(cbind(ind,vll1,vll2))

p1 <- ggplot(mat_liks1)+geom_point(data=mat_liks1,aes(ind,vll1),size=3)+ylim(-5500,-2500)+geom_hline(yintercept=-2880,linetype="dashed",colour="grey")+geom_point(data=mat_liks1,aes(ind,vll2.1),colour="darkgrey",size=1.6)+xlab("Number of States")+ylab("Mean (Median) Log-Likelihood")+ggtitle("A")
p1


# looks like the optimal number of states is 3 (where the plateau starts), so examine the results:
#   the summary contains all the estimates, and the approximate kullback-leibler divergence
summary(res_n3)
#   graphical diagnostics and confusion matrix
plot(res_n3, simulation=TRUE,true_means1,true_sigma1,true_T1,simdata1full$states)


