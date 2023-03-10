% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_poisson.R
\name{hmm_mcmc_pois}
\alias{hmm_mcmc_pois}
\title{MCMC simulation of a Hidden Markov Normal Model for Poisson model}
\usage{
hmm_mcmc_pois(
  data,
  prior_T,
  prior_betas,
  prior_alpha = 1,
  iter = 1500,
  warmup = floor(iter/1.5),
  thin = 1,
  seed = sample.int(.Machine$integer.max, 1),
  init_T = NULL,
  init_betas = NULL,
  init_alpha = NULL,
  print_params = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{data}{(numeric) data}

\item{prior_T}{(matrix) prior transition matrix}

\item{prior_betas}{(numeric) prior beta parameters}

\item{prior_alpha}{(numeric) a single prior alpha parameter. By default, \code{prior_alpha=1}}

\item{iter}{(integer) number of MCMC iterations}

\item{warmup}{(integer) number of warmup iterations}

\item{thin}{(integer) thinning parameter. By default, \code{1}}

\item{seed}{(integer) seed parameter}

\item{init_T}{(matrix) \code{optional parameter}; initial transition matrix}

\item{init_betas}{(numeric) \code{optional parameter}; initial beta parameters}

\item{init_alpha}{(numeric) \code{optional parameter}; initial alpha parameter}

\item{print_params}{(logical) \code{optional parameter}; print estimated parameters every iteration. By default, \code{TRUE}}

\item{verbose}{(logical) \code{optional parameter}; print additional messages. By default, \code{TRUE}}
}
\value{
List with following elements:
\itemize{
  \item data: data used for simulation
  \item estimates: list with various estimates
  \item idx: indices with iterations after the warmup period
  \item priors: prior parameters
  \item inits: initial parameters
  \item last_iter: list with samples from the last MCMC iteration
  \item info: list with various meta information about the object
}
}
\description{
MCMC simulation of a Hidden Markov Normal Model for Poisson model
}
\details{
Here details
}
\examples{
# Simulate Poisson-Gamma data
N <- 2^10
true_T <- rbind(c(0.95, 0.05, 0),
                c(0.025, 0.95, 0.025),
                c(0.0, 0.05, 0.95))

true_betas <- c(2, 1, 0.1)
true_alpha <- 1

simdata_full <- hmm_simulate_poisgamma_data(L = N,
                                            mat_T = true_T,
                                            betas = true_betas,
                                            alpha = true_alpha)
simdata <- simdata_full$data
hist(simdata, breaks = 40, probability = TRUE,  
     main = "Distribution of the simulated Poisson-Gamma data")
lines(density(simdata), col = "red")

# Set numbers of states to be inferred
n_states_inferred <- 3

# Set priors
prior_T <- generate_random_T(n_states_inferred)
prior_betas <- c(1, 0.5, 0.1)
prior_alpha <- 3

# Simmulation settings
iter <- 50
warmup <- floor(iter / 5) # 20 percent
thin <- 1
seed <- sample.int(10000, 1)
print_params <- FALSE # if TRUE then parameters are printed in each iteration
verbose <- FALSE # if TRUE then the state of the simulation is printed

# Run MCMC sampler
res <- hmm_mcmc_pois(data = simdata,
                     prior_T = prior_T,
                     prior_betas = prior_betas,
                     prior_alpha = prior_alpha,
                     iter = iter,
                     warmup = warmup,  
                     thin = thin,
                     seed = seed,
                     print_params = print_params,
                     verbose = verbose)
res

summary(res)# summary output can be also assigned to a variable

coef(res) # extract model estimates

# plot(res) # MCMC diagnostics
}
\references{
Here references
}
