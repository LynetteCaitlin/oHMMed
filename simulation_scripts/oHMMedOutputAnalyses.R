
# START TUTORIAL
#######################

# ********************************************************************************************
# This script produces comparative plots for...
# 1) summarizing the output of multiple runs of the MCMC algorithm 
#    on the same genomic feature with different numbers of states
# 2)+ 3) analysing the output of MCMC algorithms with fixed numbers of states
#    run on two different genomic features. In particular, this enables
#    analysis of the correlation between these features. 
# The same style of plots are presented in the paper 
#   that accompanies the first release of the package. 
# ********************************************************************************************

#--------------------------------
#  load the following libraries
#--------------------------------
 library(ggplot2)
 library(mistr)
 library(gridExtra)
 library(vcd)
 library(zoo)

# if not installed, use install(...) and then library(...)

#------------------------------
# Note: 
#-----------------------------
# The ohmmed assumption of autocorrelation can be checked via the following code:
# Set: dat <- 'ta sequence of raw data'
L=length(dat)
#     create an index for permutating the raw data:
randseq=sample(L-1,L-1)

#     run an F-test for differences in deviation from the mean 
#         in consecutive values in the sequences
#   -> THIS MUST BE SIGNIFICANT!
var.test(dat[1:(L-1)]-dat[2:L],dat[randseq]-dat[2:L])

# view the individual variances of raw and permutated sequences:
var(dat[1:(L-1)]-dat[2:L],na.rm=TRUE)
var(dat[randseq]-dat[2:L],na.rm=TRUE)



#------------------------------
# PLOT 1) 
#-----------------------------

# Assume:
#    res_n2 is the output of the class hmm.mcmc_normal or hmm.mcmc_gamma_poisson with 2 states inferred,
#    res_n3 is the output of the class hmm.mcmc_normal or hmm.mcmc_gamma_poisson with 3 states inferred,
#    etc.,... (continue for as many numbers of states as wished)

# Then: 
# create vectors of mean and median log-likelihoods from these outputs
# vll1 <- c(mean(res_n2$estimates$log_likelihood),
#           mean(res_n3$estimates$log_likelihood),
#           mean(res_n4$estimates$log_likelihood),
#           mean(res_n5$estimates$log_likelihood))
# 
# vll2 <- c(median(res_n2$estimates$log_likelihood),
#           median(res_n3$estimates$log_likelihood),
#           median(res_n4$estimates$log_likelihood),
#           median(res_n5$estimates$log_likelihood))
# 
# note: check these vector are of the same length! 
# 
# --- Panel1 - NORMAL and GAMMA-POISSON EMISSIONS
#
# Now: run the following code
ind_state <- c(res_n2$info$n_states,
               res_n3$info$n_states,
               res_n4$info$n_states,
               res_n5$info$n_states)
df_liks <- data.frame(cbind(ind_state, vll1, vll2))

# The object 'p1' stores the plot of the mean and median log-likelihoods
# for the models with increasing numbers of states
p1 <- ggplot(df_liks) + 
  geom_point(data = df_liks, aes(ind_state, vll1), size = 3) + 
  ylim(min(c(vll1, vll2)), max(c(vll1, vll2))) +
  geom_point(data = df_liks, aes(ind_state, vll2), colour = "darkgrey", size = 1.6) +
  xlab("Number of States") + 
  ylab("Mean (Median) Log-Likelihood") + 
  ggtitle("A")
p1

# --- Panel2 - NORMAL and GAMMA-POISSON EMISSIONS
# Now: select the model with appropriate number of states 'opt' - we will call it 'res_n_opt' here
# opt <- 
# res_n_opt <- 
# To store the plot of posterior means for each state, run:
nam <-  rep(paste0("state", 1:opt), each = (res_n_opt$info$iter - res_n_opt$info$warmup))
val <- as.vector(res_n_opt$samples$mean[-seq_len(res_n_opt$info$warmup), ]) # remove "warmup samples"
postm <- data.frame(name = nam, values = val)
p2 <- ggplot(postm, aes(x = name, y = values)) + 
  geom_jitter(color = "grey", size = 1.5) + 
  geom_boxplot() + 
  ggtitle("B") + 
  xlab("States") + 
  ylab("Average Value")
p2

#---- Panel3 - NORMAL EMISSIONS:
# Then: store the plot that compares the observed distribution to the theoretical 
#       distribution simulated from the estimated parameters (i.e. density plots).
#       The means of each state as well as the 95% confidence intervals 
#       are also shown in vertical lines (lines and dotted lines resp.)
state_tab <- table(res_n_opt$estimates$posterior_states)

sim_output <- unlist(lapply(1:length(state_tab), function(i) {
  stats::rnorm(state_tab[i], res_n_opt$estimates$means[i], res_n_opt$estimates$sd)
}))

dens_df <- as.data.frame(cbind(c(rep("inferred", length(sim_output)),
                                 rep("observed", length(res_n_opt$data))),
                               c(sim_output,res_n_opt$data)))

names(dens_df) <- c("data_type", "value")

p3  <- ggplot(dens_df, aes(x = as.numeric(value), color = data_type)) +
  geom_density() + 
  scale_fill_grey() +
  scale_color_manual(values = c("darkgrey", "black")) +
  labs(title = "Model Fit", x = "Values", y = "Observed (Inferred) Density") +
  geom_vline(xintercept = res_n_opt$estimates$means, color = "darkgrey", size = 0.6) +
  geom_vline(xintercept = c(res_n_opt$estimates$means) + res_n_opt$estimates$sd,
                      linetype = "dotted", color = "darkgrey", size = 0.6) +
  geom_vline(xintercept = c(res_n_opt$estimates$means) - res_n_opt$estimates$sd,
                      linetype = "dotted",  color = "darkgrey", size = 0.6) +
  guides(fill = guide_legend(title = "Density Type")) +
  ggtitle("C") + 
  theme(legend.position = "none")
p3

# ---- Panel4 - NORMAL EMISSIONS:
# Next: store the qqplot of theoretical quantiles vs quantiles simulated
#       from the mixture of posterior distributions
dims <- length(res_n_opt$estimates$means)
dist <- rep("norm", dims)  
params <- list()
for (j in 1:dims) {
  params[[j]] <- c(as.numeric(res_n_opt$estimates$means[j]), res_n_opt$estimates$sd)
}
tab_post <- table(res_n_opt$estimates$posterior_states)
weight <- as.numeric(tab_post / sum(tab_post))
resulting_mixture_of_normals <- mistr::mixdist(dist, params, weights = weight)

p4 <- mistr::QQplotgg(res_n_opt$data, resulting_mixture_of_normals, 
                      col = "black", line_col = "grey",alpha = 1) +
  ggplot2::labs(x = "Sample Quantiles", y = "Theoretical Quantiles") +
  ggplot2::theme_set(theme_bw()) + 
  ggtitle("D")
p4


#---- Panel3 - GAMMA-POISSON EMISSIONS:
# Then: store the plot that shows the observed data as a histogram
#       vertical lines demark the means of each state
data <- simdata_gp

dens_df <- as.data.frame(cbind(rep("observed", length(data)), data))
names(dens_df) <- c("data_type", "value")

p33 <- ggplot(dens_df, aes(x = as.numeric(value))) +
  geom_histogram(bins = floor(dim(table(factor(data, levels = 0:max(data))))), 
                 fill = 'grey', position = "identity") +
  scale_color_manual(values = c("black")) +
  geom_vline(xintercept = res_n_opt$estimates$means[c(1, 2)], 
             color = "black", size = 0.2) +
  labs(title = "C", x = "Number of Occurences", y = "Frequency") 
p33

#---- Panel4 - GAMMA-POISSON EMISSIONS:
# Next: simulate a theoretical distribution from the estimated parameters
#       this is shown as a fitted line over a shifted bar plot of the observed 
#       values (shifted since the bars are moved to hug the fitted line 
#       so that deviation from the theoretical distribution must be assessed
#       by checking where the base of the bars lie compared to the x-axis)

res_n_opt <- res_gp_n3
sim_output1 <- NA   
for (j in 1:500) {
  sim_output <- unlist(lapply(1:res_n_opt$inf$n_states, function(i) {
    stats::rnbinom(state_tab[i],
                   size = res_n_opt$estimates$alpha,
                   prob = res_n_opt$estimates$betas[i] / (1 + res_n_opt$estimates$betas[i]))
  }))
  sim_output1 <- c(sim_output1, sim_output)
}
sim_output1 <- sim_output1[2:length(sim_output1)]
dens_data1 <- table(factor(data, levels = 0:max(c(data, sim_output1)))) / sum(table(factor(data, levels = 0:max(c(data, sim_output1))))) 
dens_sim1 <- table(factor(sim_output1, levels = 0:max(c(data, sim_output1)))) / sum(table(factor(sim_output1, levels = 0:max(c(data, sim_output1)))))


p44 <- rootogram(as.numeric(dens_data1), 
                 as.numeric(dens_sim1), 
                 lines_gp = gpar(col = "black", lwd = 1), 
                 main = "",
                 points_gp = gpar(col = "black"), pch = "",
                 ylab = "Gene Frequencies (sqrt)",
                 xlab = "Number of Occurrences")


# The code below draws the previously stored plots in a grid in the R plotting window
# The plot can be saved/exported with additional R code

# for the normal emissions:
grid.arrange(p1, p2, p3, p4, nrow = 2)

# for the poisson emissions:
grid.arrange(p1, p2, p3, nrow = 2)
p4

#------------------------------
# PLOT 2)
#-----------------------------

# Assume:  The objects 'res_feature1' and 'res_feature2' of the class hmm.mcmc_normal or hmm.mcmc_gamma_poisson.
#   Further, we will assume in the below that the input data for both of these mcmc runs is of the same length
#   and that the entries in the input data correspond to the same windows along the genome. 
# If this is not the case, the below code must be modified (in particular for the final plot) so that each window 
#   can be assigned an inferred state according to both genomic features. 

# The first of the plots below show the input data points (which are values per genomic window) and inferred posterior means for both features along the genome. For visibility,
# one can choose to show only a fraction of the genome, for eg. chromosomes 1-6 of a total of 22. Individual genomic segments
# within the total length plotted, for eg. each chromosome in this example, can be highlighted by a bar of alternating colours
# displayed below the length of plotted genome.

# So, define the total length of the fraction of the genome to be plotted (Len),
# and then pass a vector of the window indexes contained within this:
#   Len <- 
#   windows <- 

post_means_feature1 = vector(length = Len)
for (l in 1:Len) {
  post_means_feature1[l] = sum(res_feature1$estimates$means * res_feature1$estimates$posterior_states_prob[l, ])
}

post_means_feature2 = vector(length = Len)
for (l in 1:Len) {
  post_means_feature2[l] = sum(res_feature2$estimates$means * res_feature2$estimates$posterior_states_prob[l, ])
}

# The rolling correlation the posterior means can be calculated to assess the strength of monotone correlations (Spearman correlation) :
#   pick the (overlapping) windows for which the correlation is determined:
#     width <- 
CorrDat <- as.data.frame(cbind(post_means_feature1,post_means_feature2))
rollcorr <- rollapply(CorrDat, width=width, function(x) cor(x[,1],x[,2],method="spearman"),by.column=FALSE)



# Set the limits of the y-axis to allow for space below the plotted data points for the bars to denote genomic segments
#   ylims1 <- 
#   ylims2 <- 
# Then choose a point on the y-axis within this designated space 
#   pos <- 
# Next, set the indexes of windows per each genomic segment and the length of these segments
# (such as index_a <- 1:10000, and length_a <- length(1:10000))
#   index_a <- 
#   index_b <-
#   etc.,...
#   length1 <- 
#   length2  <- 
#   etc.,...

# Next, pick colours to represent the different states for both features
#     and fill them in between the "" below (and add or remove some as necessary)
# cols1 <- c("","","","")[res_feature1$estimates$posterior_states[1:Len]]
# cols2 <- c("","","")[res_feature2$estimates$posterior_states[1:Len]]

#Then create the following grid of plots:
par(mfrow=c(3,1))
# First, the genomic landscape of the first feature:
plot(1:Len,
     res_feature1$data[1:Len],
     col = cols1,
     main = "A",
     ylab = "Feature 1 - Posterior Means",
     xaxt = "n",
     xlab = "Windows",
     ylim = ylims1)
lines(1:Len, post_means_feature1)
# The following lines can be used to designate genomic segments in alternating colours
points(cbind(index_a, rep(pos, length(index_a))), pch = "|", col = "darkgray")
points(cbind(index_b, rep(pos, length(index_b))), pch = "|", col = "lightgray")
# etc.,...

# Then, the genomic landscape of the second feature:
plot(1:Len,
  res_feature2$data[2:Len],
  col = cols2,
  main = "A",
  ylab = "Feature 2 - Posterior Means",
  xaxt = "n",
  xlab = "Windows",
  ylim = ylims2)
lines(1:Len, post_means_feature2)
# The following lines can be used to designate genomic segments in alternating colours
points(cbind(index_a, rep(pos, length(index_a))), pch = "|", col = "darkgray" )
points(cbind(index_b, rep(pos, length(index_b))), pch = "|", col = "lightgray")
# etc.,...

# Finally, the rolling correlation between the two features can be plotted:

plot(1:Len,
     rollcorr[1:Len],
     type="l",
     main="C",
     ylab="Rolling Correlation",
     xaxt="n",xlab="Windows",
     ylim=c(-1.5,1.2))

abline(h=0,col="red")
abline(h=1,col="red")
abline(h=-1,col="red")

# The following lines can be used to designate genomic segments in alternating colours
points(cbind(index_a, rep(-1.3, length(index_a))), pch = "|", col = "darkgray" )
points(cbind(index_b, rep(-1.3, length(index_b))), pch = "|", col = "lightgray")
# etc.,...

par(mfrow=c(1,1))

#------------------------------
# PLOT 3)
#-----------------------------

# The final plot cross-tabulates the assignments of each window to the two features
# This reveals the overall association between the features
Dat <- as.data.frame(cbind(res_feature1$estimates$posterior_states, res_feature2$estimates$posterior_states))
names(Dat) <- c("feature1", "feature2")

TA <- table(Dat$feature1, Dat$feature2)
barplot(TA, 
        beside = T, 
        ylab = "counts feature 1 ", 
        xlab = "states feature 2", 
        main = "")


#######################
# END TUTORIAL




