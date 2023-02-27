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
# https://github.com/LynetteCaitlin/oHMMed/UsageRecommendations.pdf
# This code is for the example simulations in the above. 
####################################################################

##### Set general MCMC parameters:
iter <- 2000               # set number of iterations; note that the default is 1000
warmup <- floor(iter / 2.5) # length of burnin is 40% of iter; note that the default is 20%
print_params <- FALSE     # parameters after each iteration will NOT be printed on the screen
verbose <- TRUE           # progress bar will be shown, as well as messages

###################################
# Appendix Example 
###################################

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
simdata1full <- hmm_simulate_poisgamma_data(L = L1,
                                            mat_T = true_T1,
                                            betas = true_betas1,
                                            alpha = true_alpha1)
# then extract the simulated data/observed sequence...:
simdata1 <- simdata1full$data
# ... and have a quick peek at the overall emission density:
hist(simdata1,breaks=50,main="")


##### Set up inference framework!:

# Set priors and initial values!

# RECOMMENDED PROCEDURE:
#   initial values and prior values should be set to the same value
#     however, for demonstration purposes we will not set initial values here
#       by default, these are then drawn from the priors
#   however, we will set all the priors according to our recommendations
#     note that these recommendations apply to overall emission densities that resemble overdispersed poisson densities
#         and user discretion is advised for overall emission densities that already resemble mixtures of bell curves

#  all prior betas will likely lie below the 'empirical overall beta':
(sum(simdata1)/length(simdata1))/((mean(simdata1)^2)/((var(simdata1)-mean(simdata1))))
#  But note that setting them near or lower than the following being a good bet:
mean(simdata1)

#  set up the priors for all runs with these recommended betas

n2_states_inferred <- 2                            # number of states to be inferred
prior2_T <- generate_random_T(n2_states_inferred) # prior transition matrix, randomly generated
prior_alpha2 <- (mean(simdata1)^2)/((var(simdata1)-mean(simdata1))/2) # recommended prior alpha! changes by number of inferred states!
prior_betas2 <- c(5,1)

n3_states_inferred <- 3
prior3_T <- generate_random_T(n3_states_inferred)
prior_alpha3 <- (mean(simdata1)^2)/((var(simdata1)-mean(simdata1))/3)
prior_betas3 <- c(5,3,1)

n4_states_inferred <- 4
prior4_T <- generate_random_T(n4_states_inferred)
prior_alpha4 <- (mean(simdata1)^2)/((var(simdata1)-mean(simdata1))/4)
prior_betas4 <- c(5,4,3,1)

n5_states_inferred <- 5
prior5_T <- generate_random_T(n5_states_inferred)
prior_alpha5 <- (mean(simdata1)^2)/((var(simdata1)-mean(simdata1))/5)
prior_betas5 <- c(5,4,3,2,1)

##### Series of oHMMed inference runs for increasing numbers of states:

res1_n2 <- hmm_mcmc_pois(data = simdata1,
                           prior_T = prior2_T,
                           prior_betas = prior_betas2,
                           prior_alpha = prior_alpha2,
                           iter = iter,
                           warmup = warmup,
                           print_params = print_params,
                           verbose = verbose)#

# Recall: it is recommended to also set: init_betas = prior_betas2,init_alpha = prior_alpha2,init_T = prior2_T

start_time <- Sys.time()
res1_n3 <- hmm_mcmc_pois(data = simdata1,
                         prior_T = prior3_T,
                         prior_betas = prior_betas3,
                         prior_alpha = prior_alpha3,
                         iter = iter,
                         warmup = warmup,
                         print_params = print_params,
                         verbose = verbose)
end_time <- Sys.time()

res1_n4 <- hmm_mcmc_pois(data = simdata1,
                         prior_T = prior4_T,
                         prior_betas = prior_betas4,
                         prior_alpha = prior_alpha4,
                         iter = iter,
                         warmup = warmup,
                         print_params = print_params,
                         verbose = verbose)

res1_n5 <- hmm_mcmc_pois(data = simdata1,
                         prior_T = prior5_T,
                         prior_betas = prior_betas5,
                         prior_alpha = prior_alpha5,
                         iter = iter,
                         warmup = warmup,
                         print_params = print_params,
                         verbose = verbose)



###--- Analyse Output and Prep plots for reading out as eps files for the overleaf document

###  
summary(res3_p_n3)
coef(res3_p_n3)
plot(res3_p_n3,simulation = TRUE,true_betas = true_betas3,true_alpha = true_alpha3,true_mat_T = true_T3p,true_states = simdata3full$states)

true_alpha3/true_betas3
as.numeric(coef(res3_p_n3)$alpha)/as.numeric(coef(res3_p_n3)$betas)

x <- conf_mat(N = L3, res = res3_p_n3, plot = TRUE)
plot_confusion_matrix(x$`Confusion Matrix`[[1]],add_sums = TRUE)

#####


vll1.1=c(mean(res3_p_n2$estimates$log_likelihood),
         mean(res3_p_n3$estimates$log_likelihood),
         mean(res3_p_n4$estimates$log_likelihood),
         mean(res3_p_n5$estimates$log_likelihood))

vll2.1=c(median(res3_p_n2$estimates$log_likelihood),
         median(res3_p_n3$estimates$log_likelihood),
         median(res3_p_n4$estimates$log_likelihood),
         median(res3_p_n5$estimates$log_likelihood))



ind <- c(2,3,4,5)
mat_liks1 <- data.frame(cbind(ind,vll1.1,vll2.1))

p1 <- ggplot(mat_liks1)+geom_point(data=mat_liks1,aes(ind,vll1.1),size=3)+geom_hline(yintercept=-122080,linetype="dashed",colour="grey")+geom_point(data=mat_liks1,aes(ind,vll2.1),colour="darkgrey",size=1.6)+xlab("Number of States")+ylab("Mean (Median) Log-Likelihood")

setEPS()
postscript("PoisLiksComp.eps", width = 5, height = 3.5,paper="special",pointsize=30)
grid.arrange(p1,nrow=1)
dev.off()


data <-res3_p_n3$data
x<- res3_p_n3
true_states <- simdata3full$states
true_means <- true_alpha3/true_betas3

states_df <- as.data.frame(cbind(1:length(data), data, x$estimates$posterior_states))
post_means <- numeric(length(data))

for (l in 1:length(data)) {
  post_means[l] <- sum(x$estimates$alpha/x$estimates$betas * x$estimates$posterior_states_prob[l, ])
}

states_df$post_means <- post_means
names(states_df) <- c("position", "data", "posterior_states","posterior_means")
states_df$posterior_states <- as.factor(states_df$posterior_states)
statesplot1_p3 <- ggplot2::ggplot(states_df, ggplot2::aes(x = position, y = data)) +
  ggplot2::geom_line(col = "grey",size=0.05) +
  ggplot2::geom_point(ggplot2::aes(colour = posterior_states), shape = 20, size = 1.5) +scale_colour_grey(start = 0.3, end = .8)+ theme_bw()+
  #ggplot2::geom_line(ggplot2::aes(x = position, y = posterior_means), size = 0.15) +
  ggplot2::guides(colour = ggplot2::guide_legend(title = "Post States")) +
  ggplot2::labs(x = "Position", y = "Data")+ggtitle("A")
#
  states_df2 <- as.data.frame(cbind(1:length(data), data, true_states))
  post_means <- (true_alpha3/true_betas3)[true_states]
  states_df2$post_means <- post_means
  names(states_df2) <- c("position", "data", "true_states","true_means")
  states_df2$true_states <- as.factor(states_df2$true_states)
  statesplot2_p3 <- ggplot2::ggplot(states_df2, ggplot2::aes(x = position, y = data)) +
    ggplot2::geom_line(col = "grey",size=0.05) +
    ggplot2::geom_point(ggplot2::aes(colour = true_states), shape = 20, size = 1.5) +scale_colour_grey(start = 0.3, end = .8)+ theme_bw()+
   # ggplot2::geom_line(ggplot2::aes(x = position, y = true_means), size = 0.15) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "True States")) +
    ggplot2::labs(x = "Position", y = "Data")+ggtitle("B")



  post_states <- x$estimates$posterior_states
  state_tab <- table(post_states)
  beta_est <- x$estimates$betas
  alpha_est <- x$estimates$alpha
  
  sim_output1 <- NA   
  for (j in 1:500) {
    sim_output <- unlist(lapply(1:length(state_tab), function(i) {
      stats::rnbinom(state_tab[i],
                     size = alpha_est,
                     prob = beta_est[i] / (1 + beta_est[i]))
    }))
    sim_output1 <- c(sim_output1,sim_output)
  }
  sim_output1 <- sim_output1[2:length(sim_output1)]
  dens_data1 <- table(factor(data,levels=0:max(c(data,sim_output1))))/sum(table(factor(data,levels=0:max(c(data,sim_output1))))) 
  dens_sim1 <- table(factor(sim_output1,levels=0:max(c(data,sim_output1))))/sum(table(factor(sim_output1,levels=0:max(c(data,sim_output1)))))
  
 
  dens_df <- as.data.frame(cbind(rep("observed", length(data)),
                                 c(data)))
  
  names(dens_df) <- c("data_type", "value")
  
  klplot <- ggplot2::ggplot(dens_df, ggplot2::aes(x = as.numeric(value))) +
    geom_histogram(bins=floor(dim(table(factor(data,levels=0:max(data))))), fill='grey',position='identity') +
    scale_color_manual(values=c("black"))+
    ggplot2::geom_vline(xintercept = true_alpha3/true_betas3, color = "white", size = 1.5)+
    ggplot2::geom_vline(xintercept = x$estimates$means, color = "black", size = 0.3) +
    ggplot2::labs(title = "Observed Counts and Inferred (and True) Means", x = "Number of Occurences", y = "Frequency")+ggtitle("C")
    
  
  setEPS()
  postscript("PoisSimsDiag2.eps", width = 10, height = 10,paper="special",pointsize=30)
  gr1a <- ggplot() + theme_void()
  gr1 <- plot_grid(klplot,gr1a,nrow=1,rel_widths=c(2/3,1/3))
  plot_grid(statesplot1_p3,statesplot2_p3,gr1,nrow=3)
  dev.off()
  
  
  
  setEPS()
  postscript("PoisDiagRoot.eps", width = 10, height = 8,paper="special",pointsize=18)
  rootogram(as.numeric(dens_data1), as.numeric(dens_sim1),lines_gp = gpar(col = "black", lwd = 2), main="", rect_gp = gpar(fill="white"),
            points_gp = gpar(col = "white"), pch = "",ylab="Frequency (sqrt)",xlab="Number of Occurrences")
  
  dev.off()
  
  
  #### change traceplot code
  
  traceplot <- function (D, family = NA, original_burnin = TRUE, original_thin = TRUE, 
                         simplify = NULL, hpd = FALSE, greek = FALSE) 
  {
    if (!is.na(family)) {
      D <- get_family(D, family = family)
    }
    if (!is.null(simplify)) {
      if (simplify > 0 & simplify < 1) {
        aD <- attributes(D)
        D <- dplyr::sample_frac(D, simplify)
        attr(D, "nChains") <- aD$nChains
        attr(D, "nParameters") <- aD$nParameters
        attr(D, "nIterations") <- aD$nIterations
        attr(D, "nBurnin") <- aD$nBurnin
        attr(D, "nThin") <- aD$nThin
        attr(D, "description") <- aD$description
      }
      else {
        stop("It is not possible to guess the simplification percentage to apply.")
      }
    }
    if (attributes(D)$nChains <= 1) {
      f <- ggplot(D, aes(x = Iteration, y = value))
    }
    else {
      f <- ggplot(D, aes(x = Iteration, y = value, colour = as.factor(Chain)))
    }
    f <- f + geom_line(alpha = 1) + scale_colour_discrete(name = "Chain")
    if (hpd) {
      ciD <- ci(D)
      min.iteration <- min(D$Iteration)
      max.iteration <- max(D$Iteration)
      offset <- min.iteration - ((max.iteration - min.iteration) * 
                                   0.01)
      f <- f + geom_segment(data = ciD, size = 0.5, inherit.aes = FALSE, 
                            aes(x = offset, xend = offset, y = low, yend = high)) + 
        geom_segment(data = ciD, size = 2, inherit.aes = FALSE, 
                     aes(x = offset, xend = offset, y = Low, yend = High))
    }
    if (!greek) {
      f <- f + facet_wrap(~Parameter, ncol = 1, scales = "free")
    }
    else {
      f <- f + facet_wrap(~Parameter, ncol = 1, scales = "free", 
                          labeller = label_parsed)
    }
    t_format <- function(x) {
      return(((x - 1) * attributes(D)$nThin) + attributes(D)$nThin)
    }
    b_format <- function(x) {
      return(x + attributes(D)$nBurnin)
    }
    bt_format <- function(x) {
      return(attributes(D)$nBurnin + (((x - 1) * attributes(D)$nThin) + 
                                        attributes(D)$nThin))
    }
    if (original_burnin & !original_thin) {
      f <- f + scale_x_continuous(labels = b_format)
    }
    else if (!original_burnin & original_thin) {
      f <- f + scale_x_continuous(labels = t_format)
    }
    else if (original_burnin & original_thin) {
      f <- f + scale_x_continuous(labels = bt_format)
    }
    else {
      f <- f
    }
    return(f)
  }
  
  #### end change traceplot code
  
x <- res3_p_n3  
data <- x$data
idx <- x$idx
  
  
# Diagnostics betas
all_betas <- convert_to_ggmcmc(x, pattern = "beta")
n_betas <- attributes(all_betas)$nParameters
facet_means <- ggplot2::facet_wrap(~ Parameter, ncol = floor(n_betas / 2), scales = "free")
mdens <- ggmcmc::ggs_density(all_betas) + facet_means
mtrace <- traceplot(all_betas) + facet_means


# Diagnostics transitions
all_T <- convert_to_ggmcmc(x, pattern = "T")
n_t <- attributes(all_T)$nParameter
facet_t <- ggplot2::facet_wrap(~ Parameter, ncol = sqrt(n_t), scales = "free")
labels_t <- ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = "."))
Ttrace <- traceplot(all_T) + facet_t + labels_t + ggplot2::labs(x = "Iteration", y = "Value")
Tdens <- ggmcmc::ggs_density(all_T) + facet_t + labels_t + ggplot2::labs(x = "Value", y = "Density")

# Diagnostics alpha

  df_alpha <- convert_to_ggmcmc(x, "alpha")
  sdtrace <- traceplot(df_alpha) + ggplot2::labs(x = "Iteration", y = "Value")
  sddens <- ggmcmc::ggs_density(df_alpha) + ggplot2::labs(x = "Value", y = "Density")

#Diagnostics mean
all_means <- convert_to_ggmcmc(x, pattern = "means")
n_means <- attributes(all_means)$nParameters
facet_me <- ggplot2::facet_wrap(~ Parameter, ncol = floor(n_means / 2), scales = "free")
medens <- ggmcmc::ggs_density(all_means) + facet_me
metrace <- traceplot(all_means) + facet_me

# Likelihood trace
lltrace <- x$estimates$log_likelihood[idx]
lltrace_df <- as.data.frame(cbind(idx,lltrace))
names(lltrace_df) <- c("iteration", "log_likelihood")

llplot <- ggplot2::ggplot(lltrace_df, ggplot2::aes(x = iteration, y = log_likelihood)) +
  ggplot2::geom_line() +
  ggplot2::labs(x = "Iteration", y = "Log-likelihood")


setEPS()
postscript("PoisSimsDiag1a.eps", width = 12, height = 4,paper="special",pointsize=30)
plot_grid(llplot,mdens,mtrace,sddens,sdtrace,nrow=1, rel_widths = c(1/3,0.5/3,0.5/3,0.5/3,0.5/3))
dev.off()

setEPS()
postscript("PoisSimsDiag1b.eps", width = 10, height = 8,paper="special",pointsize=30)
plot_grid(Tdens,medens,Ttrace,metrace,nrow=2, rel_widths = c(2/3,1/3))
dev.off()

######################################### !!!

###################################
# Appendix Checks 
###################################

L4 <-
#(at least 2^11 but better over 2^13)
  
true_T4p <- rbind(c(0.99, 0.01, 0),
                  c(0.01, 0.98, 0.01),
                  c(0.0, 0.01, 0.99))

true_betas4 <- 1 / (c(0.2, 1.5, 9))
true_alpha4 <- 1


simdata4full <- hmm_simulate_poisgamma_data(L = L4,
                                            mat_T = true_T4p,
                                            betas = true_betas4,
                                            alpha = true_alpha4)
simdata4 <- simdata4full$data

start_time <- Sys.time()
res4_p_n3 <- hmm_mcmc_pois(data = simdata4,
                           prior_T = true_T4p,
                           prior_betas = true_betas4,
                           prior_alpha = true_alpha4,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

end_time <- Sys.time()

# iter 1000
# Time difference of 1.67714 mins 2^12 (4096)-100.6284
# Time difference of 3.364693 mins 2^13 (8192)- 201.88158
#Time difference of 6.688331 mins for 2^14 (16384)- 401.29986
#Time difference of 13.46753 mins for 2^15 (32768) - 808.0518  
# Time difference of 26.98619 mins - 1619.1714

# compare: iter  (time/5000)*10

setEPS()
postscript("PoisSysTime.eps", width = 6, height = 4,paper="special",pointsize=10)
plot(c(2^12,2^13,2^14,2^15,2^16),c((100.6284/1000)*10,(201.88158/1000)*10,(401.29986/1000)*10,(808.0518/1000)*10,(1619.1714/1000)*10),xlab="number of windows", ylab="system time per 10 iter (seconds)",pch="x",type="b",xlim=c(0,66000),ylim=c(0,17))
dev.off()


#####
output_pois_sims<-list(res3_p_n2,res3_p_n3,res3_p_n4,res3_p_n5)
save(output_pois_sims,file="Output_pois_sims.Rdata")
  
