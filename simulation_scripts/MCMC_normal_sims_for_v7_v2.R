setwd("~/Documents/Documents - Lynette’s MacBook Pro/Research/Active/MalariaProj/optimized_hmm_scripts")

#### load code
####  as well as dependencies

library(ggmcmc)
library(ggplot2)
library(gridExtra)
library(cvms)
library(gtools)
library(mistr)
library(cowplot)

source("MCMC_normal_v7.R")

#### general MCMC parameters for all the simulations
iter <- 600
warmup <- floor(iter / 5) # 20%
thin <- 1
print_params <- FALSE
verbose <- TRUE

###################################
# Appendix Example 
###################################

L0=2^13
sigma1 <- 0.195
means1 <- c(0.55,1,1.9)
true_T1 <- rbind(c(0.7, 0.3, 0),
                 c(0.3, 0.6, 0.1),
                 c(0.0, 0.35, 0.65))
simdata1full <- hmm_simulate_normal_data(L0,true_T1,means1,sigma1)
simdata1 <- simdata1full$data

# -> quick view
plot(density(simdata1), main="")


#  Set different priors: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----

prior_sm_sd <- 0.02

n2_states_inferred <- 2
prior2_T <- generate_random_T(n2_states_inferred)
prior2_means <- c(0,3)
prior2_sd <- 2.5/2

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

#  Run different MCMCs: ----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----—----

res1_n2 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior2_T,
                           prior_means = prior2_means,
                           prior_sd = prior2_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

res1_n3 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior3_T,
                           prior_means = prior3_means,
                           prior_sd = prior3_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

res1_n4 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior4_T,
                           prior_means = prior4_means,
                           prior_sd = prior4_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

res1_n5 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior5_T,
                           prior_means = prior5_means,
                           prior_sd = prior5_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

res2_n2 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior2_T,
                           prior_means = prior2_means,
                           prior_sd = prior_sm_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

res2_n3 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior3_T,
                           prior_means = prior3_means,
                           prior_sd = prior_sm_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

res2_n4 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior4_T,
                           prior_means = prior4_means,
                           prior_sd = prior_sm_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)

res2_n5 <- hmm_mcmc_normal(data = simdata1,
                           prior_T = prior5_T,
                           prior_means = prior5_means,
                           prior_sd = prior_sm_sd,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)


###--- Analyse Output and Prep plots for reading out as eps files for the overleaf document

# check all summaries.... -> res1_n3 is the 'correct' one (as is technically res2_n4)

x <- conf_mat(N=L0,res=res1_n3,plot=TRUE) # SHOULD WORK BUT DOES NOT???? SO...
plot_confusion_matrix(x$`Confusion Matrix`[[1]],add_sums = TRUE)



vll1.1=c(mean(res1_n2$estimates$log_likelihood),
       mean(res1_n3$estimates$log_likelihood),
       mean(res1_n4$estimates$log_likelihood),
       mean(res1_n5$estimates$log_likelihood))

vll2.1=c(median(res1_n2$estimates$log_likelihood),
       median(res1_n3$estimates$log_likelihood),
       median(res1_n4$estimates$log_likelihood),
       median(res1_n5$estimates$log_likelihood))

vll1.2=c(mean(res2_n2$estimates$log_likelihood),
         mean(res2_n3$estimates$log_likelihood),
         mean(res2_n4$estimates$log_likelihood),
         mean(res2_n5$estimates$log_likelihood))

vll2.2=c(median(res2_n2$estimates$log_likelihood),
         median(res2_n3$estimates$log_likelihood),
         median(res2_n4$estimates$log_likelihood),
         median(res2_n5$estimates$log_likelihood))




ind <- c(2,3,4,5)
mat_liks1 <- data.frame(cbind(ind,vll1.1,vll2.1))
mat_liks2 <- data.frame(cbind(ind,vll1.2,vll2.2))

p1 <- ggplot(mat_liks1)+geom_point(data=mat_liks1,aes(ind,vll1.1),size=3)+ylim(-5500,-2500)+geom_hline(yintercept=-2880,linetype="dashed",colour="grey")+geom_point(data=mat_liks1,aes(ind,vll2.1),colour="darkgrey",size=1.6)+xlab("Number of States")+ylab("Mean (Median) Log-Likelihood")+ggtitle("A")
p2 <- ggplot(mat_liks2)+geom_point(data=mat_liks2,aes(ind,vll1.2),size=3)+ylim(-5500,-2500)+geom_hline(yintercept=-2895,linetype="dashed",colour="grey")+geom_point(data=mat_liks1,aes(ind,vll2.2),colour="darkgrey",size=1.6)+xlab("Number of States")+ylab("")+ggtitle("B")


setEPS()
postscript("NormalLiksComp.eps", width = 10, height = 7,paper="special",pointsize=30)
grid.arrange(p1,p2,nrow=1)
dev.off()



data <-res1_n3$data
x<- res1_n3
true_states <- simdata1full$states
true_means <- means1

states_df <- as.data.frame(cbind(1:length(data), data, x$estimates$posterior_states))
post_means <- numeric(length(data))

for (l in 1:length(data)) {
  post_means[l] <- sum(x$estimates$means * x$estimates$posterior_states_prob[l, ])
}

states_df$post_means <- post_means
names(states_df) <- c("position", "data", "posterior_states", "posterior_means")
states_df$posterior_states <- as.factor(states_df$posterior_states)
statesplot_res1_3 <- ggplot2::ggplot(states_df, ggplot2::aes(x = position, y = data)) +
  ggplot2::geom_line(col = "grey",size=0.05) +
  ggplot2::geom_point(ggplot2::aes(colour = posterior_states), shape = 20, size = 1.5) +scale_colour_grey(start = 0.3, end = .8)+ theme_bw()+
  #ggplot2::geom_line(ggplot2::aes(x = position, y = posterior_means), size = 0.01) +
  ggplot2::guides(colour = ggplot2::guide_legend(title = "Post State")) +
  ggplot2::labs(x = "", y = "Data")+ggtitle("A")

states_df2 <- as.data.frame(cbind(1:length(data), data, true_states))
post_means <- true_means[true_states]
states_df2$post_means <- post_means
names(states_df2) <- c("position", "data", "true_states", "posterior_means")
states_df2$true_states <- as.factor(states_df2$true_states)
statesplot2_res1_3 <- ggplot2::ggplot(states_df2, ggplot2::aes(x = position, y = data)) +
  ggplot2::geom_line(col = "grey",size=0.05) +
  ggplot2::geom_point(ggplot2::aes(colour = true_states), shape = 20, size = 1.5) +scale_colour_grey(start = 0.3, end = .8)+ theme_bw()+
  #ggplot2::geom_line(ggplot2::aes(x = position, y = posterior_means), size = 0.01) +
  ggplot2::guides(colour = ggplot2::guide_legend(title = "True State")) +
  ggplot2::labs(x = "Position", y = "Data")+ggtitle("C")

post_states <- x$estimates$posterior_states
state_tab <- table(post_states)
estim_means <- x$estimates$means
estim_sd <- x$estimates$sd

sim_output <- unlist(lapply(1:length(state_tab), function(i) {
  stats::rnorm(state_tab[i], estim_means[i], estim_sd)
}))
dens_df <- as.data.frame(cbind(c(rep("inferred", length(sim_output)),
                                 rep("observed", length(data))),
                               c(sim_output, data)))

names(dens_df) <- c("data_type", "value")

dens_res1_3  <- ggplot2::ggplot(dens_df, ggplot2::aes(x = as.numeric(value),color =data_type)) +
  ggplot2::geom_density() + scale_fill_grey()+ scale_color_manual(values=c("darkgrey", "black"))+
  ggplot2::labs(title = "Model Fit", x = "", y = "Observed (Inferred) Density") +
  ggplot2::geom_vline(xintercept = x$estimates$means, color = "darkgrey", size = 0.6) +
  ggplot2::geom_vline(xintercept = c(x$estimates$means) + x$estimates$sd,
                      linetype = "dotted", color = "darkgrey", size = 0.6) +
  ggplot2::geom_vline(xintercept = c(x$estimates$means) - x$estimates$sd,
                      linetype = "dotted",  color = "darkgrey", size = 0.6) +
  ggplot2::guides(fill = ggplot2::guide_legend(title = "Density Type"))+ggtitle("B")+theme(legend.position = "none")

dims <- length(x$estimates$means)
dist <- rep("norm", dims)    #### NOTE: distribution here will have to be changed for count version!!
params <- list()
for (j in 1:dims) {
  params[[j]] <- c(as.numeric(x$estimates$means[j]), x$estimates$sd)
}
tab_post <- table(x$estimates$posterior_states)
weight <- as.numeric(tab_post / sum(tab_post))
resulting_mixture_of_normals <- mistr::mixdist(dist, params, weights = weight)

qqplot_res_1_3 <- mistr::QQplotgg(data, resulting_mixture_of_normals, col = "black", line_col = "grey",alpha=1) +
  ggplot2::labs(x = "Resulting mixture of normals", y = "Data") +
  ggplot2::theme_get()+ggtitle("D")


setEPS()
postscript("NormSimsDiag2.eps", width = 10, height = 10,paper="special",pointsize=30)
plot_grid(statesplot_res1_3,dens_res1_3,statesplot2_res1_3,qqplot_res_1_3,rel_widths = c(2.5/4, 1.5/4, 2.5/4,1.5/4))
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

x <- res1_n3
info <- x$info
data <- x$data
idx <- x$idx

# Diagnostics mean
all_means <- convert_to_ggmcmc(x, pattern = "mean")
n_means <- attributes(all_means)$nParameters
facet_means <- ggplot2::facet_wrap(~ Parameter, ncol = floor(n_means / 2), scales = "free")
mdens <- ggmcmc::ggs_density(all_means) + facet_means
mtrace <- traceplot(all_means) + facet_means

# Diagnostics transitions
all_T <- convert_to_ggmcmc(x, pattern = "T")
n_t <- attributes(all_T)$nParameter
facet_t <- ggplot2::facet_wrap(~ Parameter, ncol = sqrt(n_t), scales = "free")
labels_t <- ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = "."))
Ttrace <-traceplot(all_T) + facet_t + labels_t + ggplot2::labs(x = "Iteration", y = "Value")
Tdens <- ggmcmc::ggs_density(all_T) + facet_t + labels_t + ggplot2::labs(x = "Value", y = "Density")

# Diagnostics sd
df_sigma <- convert_to_ggmcmc(x, "sigma")
sdtrace <-traceplot(df_sigma) + ggplot2::labs(x = "Iteration", y = "Value")
sddens <- ggmcmc::ggs_density(df_sigma) + ggplot2::labs(x = "Value", y = "Density")

# Likelihood trace
lltrace <- x$estimates$log_likelihood[idx]
lltrace_df <- as.data.frame(cbind(idx,lltrace))
names(lltrace_df) <- c("iteration", "log_likelihood")

llplot <- ggplot2::ggplot(lltrace_df, ggplot2::aes(x = iteration, y = log_likelihood)) +
  ggplot2::geom_line() +
  ggplot2::labs(x = "Iteration", y = "Log-likelihood")


setEPS()
postscript("NormSimsDiag1a.eps", width = 10, height = 6,paper="special",pointsize=30)
plot_grid(llplot,mdens,mtrace,nrow=1, rel_widths = c(1/2,1/4,1/4))
dev.off()

setEPS()
postscript("NormSimsDiag1b.eps", width = 10, height = 8,paper="special",pointsize=30)
plot_grid(Tdens,sddens,Ttrace,sdtrace,nrow=2, rel_widths = c(2/3,1/3))
dev.off()


###################################
# Appendix Checks 
###################################
iter <- 600
warmup <- floor(iter / 5) # 20%
thin <- 1

#L_a = 
# set to something
##(at least 2^9 but better over 2^10)
true_T_mat <- rbind(c(1-0.01,0.01,0),
               c(0.01,1-0.02,0.01),
               c(0,0.01,1-0.01))
true_means <- c(-1, 0.0, 1)
true_sd  <- sqrt(0.1)

simdata_a_full <- hmm_simulate_normal_data(L_a,true_T_mat,true_means,true_sd)
simdata_a <- simdata_a_full$data


start_time <- Sys.time()
res_a_n3 <- hmm_mcmc_normal(data = simdata_a,
                           prior_T = true_T_mat,
                           prior_means = true_means,
                           prior_sd = true_sd/3,
                           iter = iter,
                           warmup = warmup,
                           thin = thin,
                           print_params = print_params,
                           verbose = verbose)
end_time <- Sys.time()

# iter 600
#Time difference of 57.42708 secs for 2^12 (4096)
#Time difference of 1.898766 mins for 2^13 (8192) - 113.92596 secs
#Time difference of 3.965682 mins for 2^14 (16384)- 237.94092
# Time difference of 7.717089 mins for 2^15 (32768)- 463.02534
# 
# iter 1000
#Time difference of 12.87336 mins for 2^15 (32768) - 922.7958 
# Time difference of 25.66149 mins - 1539.6894

setEPS()
postscript("OverallSysTime.eps", width = 6, height = 4,paper="special",pointsize=10)
plot(c(2^12,2^13,2^14,2^15,2^16),c((57.42708/600)*10,(113.92596/600)*10,(237.94092/600)*10,(463.02534/600)*10, (1539.6894/1000)*10),xlab="number of windows", ylab="system time per 10 iter (seconds)",pch="x",type="b",xlim=c(0,66000),ylim=c(0,17))
# add pois sys time
lines(c(2^12,2^13,2^14,2^15,2^16),c((100.6284/1000)*10,(201.88158/1000)*10,(401.29986/1000)*10,(808.0518/1000)*10,(1619.1714/1000)*10),xlab="number of windows", ylab="system time per 10 iter (seconds)",pch="x",type="b",xlim=c(0,66000),ylim=c(0,17),col="grey")
dev.off()

output_norm_sims_run1<-list(res1_n2,res1_n3,res1_n4,res1_n5)
output_norm_sims_run2<-list(res2_n2,res2_n3,res2_n4,res2_n5)
save(output_norm_sims_run1,file="Output_norm_sims_run1.Rdata")
save(output_norm_sims_run2,file="Output_norm_sims_run2.Rdata")

