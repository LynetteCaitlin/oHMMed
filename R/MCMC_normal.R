### CHANGES v6:
# gtools replaces MCMCprecision
#  - 2x faster for random number generations)
#  - only one dependency for Dirichlet distribution
# main hmm_mcmc_normal function returns original dataset. It is then used
# in summary and plot functions.



### CHANGES prior v6:
# Function names:
#   - eigenSystem --> eigen_system
#   - KullbackLeibler --> kullback_leibler
#   - sampleMeansSD --> sample_means_sd_
#   - sampleStates_normal_ --> sample_states_normal_
#   - sampleT --> sample_T_
#
# New helpers:
#   - get_pi
#   - get_pi_
#   - cap_floor_
#   - init_hmm_mcmc_normal (initializes the mcmc normal algorithm and checks inputs)
#
# Changes in style
#   - changed order of inputs in functions: optional parameters that have NULL
#     go to the end of the function call
#   - all internal functions end with "_"
#   - all functions that are available to a user
#     don't end with "_"
#   - what functions should be exported (available to the user) and which not
#     will be decided in later stage of the project (easy to adapt code).
#   - Example: "get_pi" should be available to users
#     and should be "safe" to use (sanity checks).
#     However, those checks are not necessary for internal use.
#     Therefore internally "get_pi_" is used without checks
#     and users use "get_pi" which uses "get_pi_" internally and adds sanity checks.
#
# Changes in names of inputs/outputs and variables
#   - "vmeans" and "vMeans" changed to "means" everywhere.
#   - "matT" changed to "mat_T" everywhere
#
# Other changes:
#   - changed order and names of inputs to "hmm_mcmc_normal"
#     all follow "underscore_case"
#   - new inputs in "hmm_mcmc_normal": chains, iter, warmup, thin, init, ...
#
#
## DETAILED CHANGES:
# "posterior_probabilities_normal": -------------------------------------------------------------------------
# MM: few adjustments have been made to make the function robust.
#
# Main story: if the prior values (also initial values) for means and the variance
#    are far, far away from the actual sample (central) moments then
#    the probabilities calculated as:
#        dnorm(x = data[i], mean = means, sd = sdev)
#    are going to be 0 due to floating point precision.
#    This causes a chain of adverse events leading to NAs, NaNs, etc, which
#    cause issues.
#
# More details:
#    The smallest possible number is about 5e-324 (see ?.Machine)
#    and if one subtracts a tiny bit from that then the result is going to be 0.
#    Check:
#        1 * 5e-324 # = 4.940656e-324
#        1 * 2e-324 # = 0
#        1 * 5e-325 # = 0
#    The probability density is non-zero for all real numbers, however,
#    it is common that the R function returns zero due to the precision
#    problem described above. Example:
#        dnorm(1000, mean = -1e7, sd = 0.1)
#
# Remedy:
#    To overcome this issue a floor is implemented
#    defined at 2.225074e-308 == .Machine$double.xmin.
#    The density value is then as close as possible to its "true value",
#    which is obviously lower than 2.225074e-308 but not 0.
#
#    2.225074e-308 still allows performing all arithmetic operations.
#    1 / 2.225074e-308
#
# Another potential issue: 0 * Inf = NaN
#    This problem can arise when the vector pi contains 0 values
#    and dnorm produces Inf. Example:
#       dnorm(5, mean = 5, sd = 0)
#    In order to avoid the issue dnorm is capped at .Machine$double.xmax
#    which is 1.797693e+308.
#    0 * Inf --> NaN
#    0 * 1.797693e+308 --> 0
#
# Capping and flooring probabilities/densities at values close
# to limits is a reasonable technical tricks that assures stability,
# robustness (or makes calculation possible at all)
# without influencing results too much.



# "sample_T_": ------------------------------------------------------------------------------------------------
# MM: frequently a row sum is equal to 0.9999999999999998889777,
# which is not precisely 1 and therefore 0.9999999999999998889777 == 1 --> FALSE.
# This is caused by floating point precision. See typical example:
# sqrt(2)^2 == 2 # --> FALSE
# A new helper "is_row_sum_one_" tolerates a discrepancy of
# at most tol=.Machine$double.eps^(1/2) by default and is therefore robust.
# It is used throughout the main MCMC function for sanity checks.
# --> "sample_T_()" doesn't need to be repeated
# https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal

# "sample_states_normal_": -----------------------------------------------------------------------------------
# Improvements:
# 1) states are returned as factors (factor(states, levels = 1:n_states))
#    table(states) always includes all states, even if they are not found in
#    the vector (levels are stored as meta-attributes). Example below:
#      x <- factor(1:10)
#      attributes(x)
#      table(x[1:2])
# 2) dnorm is used in a vectorized way - probabilities are calculated only once
# 3) sample() is replaced with two-three fold more efficient alternative


# TODO:
# - I believe, sample_T_ and sample_means_sd should be available to the user



# Check convergence: https://www.rdocumentation.org/packages/ConvergenceConcepts/versions/1.2.2/topics/check.convergence
# Convergence after x step
# https://stackoverflow.com/questions/38143932/how-to-find-when-a-matrix-converges-with-a-loop
#
# Convergent matrix: https://en.wikipedia.org/wiki/Convergent_matrix



###############
###############
# NEW:
# How do we check that the users have the all the dependencies
# (i.e. dependent libraries for all our functions) installed?
# Can one check for it when our library is loaded ?
###############
###############




# ******************************************************************************************************-----
# GENERAL HELPERS: ------------------------------------------------------------------------------------------
# ******************************************************************************************************-----


#' @keywords internal
get_pi_ <- function(mat_T) {
  ev <- eigen_system(mat_T)
  Re(ev$forwards[1, ]) / sum(Re(ev$forwards[1, ]))
}


#' @keywords internal
cap_floor_ <- function(x, cap = Inf, floor = -Inf) {
  x[x > cap] <- cap
  x[x < floor] <- floor
  x
}


#' @keywords internal
print_progress_ <- function(tt, N, id = NULL, verbose = TRUE) {
  
  xx <- NULL
  if (!is.null(id)) {
    xx <- paste0("Chain ", id, ": ")
  }
  if (verbose) {
    if (tt %% round(N / 10) == 0) {
      message(paste0(xx, round(floor(tt / N * 100), -1), "%"))
    }
  }
}


#' @keywords internal
is_row_sum_one_ <- function(mat = NULL, tol = .Machine$double.eps^0.5) {
  rs <- rowSums(mat)
  abs(rs - 1) < tol
}


#' @keywords internal
get_mat_T_ <- function(u, l) {
  
  DIM <- length(u) + 1
  
  if ((length(u) != length(l)) | (length(u) < 1)) {
    warning("Length of `u` and `l` must be identical, greater than 0, and one less than the dimension of the matrix.", call. = FALSE)
  }
  mat_T <- matrix(rep(0, DIM * DIM), nrow = DIM)
  mat_T[1,1] <- 1.0 - u[1]
  mat_T[1,2] <- u[1]
  mat_T[DIM,DIM - 1] <- l[DIM - 1]
  mat_T[DIM,DIM] <- 1.0 - l[DIM - 1]
  
  if (DIM > 2) {
    for (i in 2:(DIM - 1)) {
      mat_T[i,i - 1] <- l[i - 1]
      mat_T[i,i + 1] <- u[i]
      mat_T[i,i] <- 1 - l[i - 1] - u[i]
    }
  }
  mat_T
}


# #' @keywords internal
# sample_integers_ <- function(mat) {
#     u <- stats::runif(nrow(mat))
#     mat <- matrixStats::rowCumsums(mat)
#     max.col(mat - u > 0, ties.method = "first")
# }


#' Converts MCMC samples into \code{ggmcmc} format
#'
#' This helper function converts MCMC samples into \code{ggmcmc} format
#'
#' @param x (mcmc_hmm_*) MCMC HMM object
#'
#' @param pattern (character) pattern(s) with model parameters to be included in the output
#'
#' @param include_warmup (logical) include warmup samples. By default \code{FALSE}
#'
#' @return
#' data.frame compatible with functions from the \code{ggmcmc} package
#'
#' @export
#'
#' @examples
#' # TODO

convert_to_ggmcmc <- function(x,
                              pattern = c("mean", "sigma", "beta", "alpha","pois_means", "T"),
                              include_warmup = FALSE) {
  
  info <- x$info
  iter <- info$iter
  n_states <- length(x$priors[[1]])
  s <- x$samples
  
  ind_mat <- as.vector(outer(1:n_states, 1:n_states, function(x,y) paste0("T[", x, ",", y, "]")))
  mat_T <- t(sapply(1:iter, function(i) { as.numeric(s$mat_T[, ,i])}))
  colnames(mat_T) <- ind_mat
  
  res <- cbind(s$means, s$betas, "sigma[1]" = s$sd, "alpha[1]" = s$alpha, s$pois_means,mat_T)
  
  nIterations <- iter
  nBurnin <- 0
  
  res <- reshape2::melt(res)
  res$Chain <- 1
  colnames(res) <- c("Iteration", "Parameter", "value", "Chain")
  res <- res[ ,c(1,4,2,3)]
  
  if (!include_warmup) {
    res <- res[res$Iteration %in% x$idx, ]
    nIterations <- iter - info$warmup
  }
  pat <- gsub("[", replacement = "\\[", x = pattern, fixed = TRUE)
  pat <- gsub("]", replacement = "\\]", x = pat, fixed = TRUE)
  pat <- paste0(pat, collapse = "|")
  
  res <- res[grepl(pat, res$Parameter), ]
  
  n_params <- length(unique(res$Parameter))
  
  attr(res, "nChains") <- 1
  attr(res, "nParameters") <- n_params
  attr(res, "nIterations") <- nIterations
  attr(res, "nBurnin") <- nBurnin
  attr(res, "nThin") <- info$thin
  attr(res, "description") <- paste(info$model_name, info$date)
  res
}


#' Generate a random transition matrix 
#'
#' This helper function generates a transition matrix at random for testing purposes
#'
#' @param n (integer) dimension of a transition matrix
#'
#' @details
#' Uniform random numbers \eqn{[0,1]} are used to fill the matrix. Rows are then
#' normalized.
#'
#' @return
#' random \code{n x n} transition matrix
#'
#' @export
#'
#' @examples
#' mat_T <- generate_random_T(3)
#' mat_T
#' rowSums(mat_T)

generate_random_T <- function(n = 3) {
  m <- diag(stats::runif(n))
  for (i in 2:n) {
    m[i,i-1] <- stats::runif(1)
    m[i-1,i] <- stats::runif(1)
  }
  for (i in 1:n) {
    m[i, ] <- m[i, ] / sum(m[i, ])
  }
  m
}


#' Calculate a continuous approximation of the Kullback-Leibler divergence
#'
#' @param p (numeric) probabilities
#'
#' @param q (numeric) probabilities
#'
#' @details
#' The continuous approximation of the Kullback-Leibler divergence
#' is calculated as follows:
#' \deqn{
#'   \frac{1}{n}\sum_{i=1}^n\big[\log(p_i) p_i - log(q_i) p_i \big]
#' }
#'
#' @return
#' Numeric vector
#'
#' @export
#'
#' @examples
#' # Simulate n normally distributed variates
#' n <- 1000
#' dist1 <- rnorm(n)
#' dist2 <- rnorm(n, mean = 0, sd = 2)
#' dist3 <- rnorm(n, mean = 2, sd = 2)
#' 
#' # Estimate probability density functions
#' pdf1 <- density(dist1)
#' pdf2 <- density(dist2)
#' pdf3 <- density(dist3)
#' 
#' # Visualise PDFs
#' plot(pdf1, main = "PDFs", col = "red", xlim = range(dist3))
#' lines(pdf2, col = "blue")
#' lines(pdf3, col = "green")
#' 
#' # PDF 1 vs PDF 2
#' kullback_leibler_cont_appr(pdf1$y, pdf2$y)
#' 
#' # PDF 1 vs PDF 3
#' kullback_leibler_cont_appr(pdf1$y, pdf3$y)
#' 
#' # PDF 2 vs PDF 2
#' kullback_leibler_cont_appr(pdf2$y, pdf3$y)

kullback_leibler_cont_appr <- function(p, q) {
  
  if (!is.numeric(p)) {
    stop("kullback_leibler_cont_appr(): `p` must be numeric", call. = FALSE)
  }
  if (!is.numeric(q)) {
    stop("kullback_leibler_cont_appr(): `p` must be numeric", call. = FALSE)
  }
  if (length(p) != length(q)) {
    stop("kullback_leibler_cont_appr(): `p` must be the same length as `q`", call. = FALSE)
  }
  (sum(log(p) * p) - sum(log(q) * p)) / length(p)
}


#' Get the prior probability of states
#'
#' Calculate the prior probability of states that correspond to the stationary
#' distribution of the transition matrix T
#'
#' @param mat_T (matrix) transition matrix
#'
#' @details
#' It is assumed that the prior probability of states corresponds
#' to the stationary distribution of the transition matrix \eqn{T},
#' denoted with \eqn{\pi} and its entries with \eqn{\pi_i=Pr(\theta_{l-1}=i)}.
#'
#' @return
#' A numeric vector
#'
#' @export
#'
#' @examples
#' T_mat <- rbind(c(1-0.01,0.01,0),
#'                c(0.01,1-0.02,0.01),
#'                c(0,0.01,1-0.01))
#' T_mat
#' get_pi(T_mat)

get_pi <- function(mat_T = NULL) {
  
  if (is.null(mat_T)) {
    stop("get_pi(): `mat_T` is not specified", call. = FALSE)
  }
  
  if (!is.matrix(mat_T) || diff(dim(mat_T)) != 0) {
    stop("get_pi(): `mat_T` must be a square matrix", call. = FALSE)
  }
  
  if (any(rowSums(mat_T) != 1)) {
    stop("get_pi(): rows in the transition matrix `mat_T` must sum up to 1", call. = FALSE)
  }
  
  get_pi_(mat_T)
}


#' Calculate eigenvalues and eigenvectors.
#'
#' This helper function returns the eigenvalues in lambda and the left and right eigenvectors in forwards and backwards.
#'
#' @param mat (matrix) a square matrix
#'
#' @return
#' a list with three elements:
#' \itemize{
#' \item lambda: eigenvalues
#' \item forwards: left eigenvector
#' \item backwards: right eigenvector
#' }
#'
#' @export
#'
#' @examples
#' mat_T0 <- rbind(c(1-0.01,0.01,0),
#'                c(0.01,1-0.02,0.01),
#'                c(0,0.01,1-0.01))
#' eigen_system(mat_T0)

eigen_system <- function(mat) {
  
  es <- eigen(mat)
  backwards <- es$vectors
  N <- length(mat[1, ])
  
  for(i in 1:length(mat[ ,1])) {
    backwards[ ,i] <- N * backwards[ ,i] / sum(sqrt(backwards[ ,i] * backwards[ ,i]))
  }
  
  forwards <- solve(backwards)
  lambda <- es$values
  list("lambda" = lambda, "forwards" = forwards, "backwards" = backwards)
}


#' Hidden Markov Model simulation with normal data.
#'
#' Simulate a Hidden Markov Model based on normal data.
#'
#' @param L (integer) number of simulations
#'
#' @param mat_T (matrix) a square matrix with the initial state
#'
#' @param means (numeric) \code{mean} parameter in \code{\link{rnorm}} for emission probabilities
#'
#' @param sigma (numeric) \code{sd} parameter in \code{\link{rnorm}} for emission probabilities
#'
#' @return
#' returns a data vector "data", the "true" hidden states "states" used to generate the data vector
#' and prior probability of states "pi".
#'
#' @export
#'
#' @examples
#' mat_T0 <- rbind(c(1-0.01,0.01,0),
#'                c(0.01,1-0.02,0.01),
#'                c(0,0.01,1-0.01))
#' L <- 2^5
#' means0 <- c(-1,0,1)
#' sigma0 <- 1
#'
#' sim_data <- hmm_simulate_normal_data(L = L, mat_T = mat_T0, means = means0, sigma = sigma0)
#' plot(density(sim_data$data), main = "Density")
#' sim_data

hmm_simulate_normal_data <- function(L, mat_T, means, sigma) {
  
  if (length(sigma) != 1) {
    stop("hmm_simulate_normal_data(): standard deviation `sigma` must be of length 1", call. = FALSE)
  }
  
  if (!is.matrix(mat_T)) {
    stop("hmm_simulate_normal_data(): `mat_T` must be a numeric matrix", call. = FALSE)
  }
  
  if (any(is_row_sum_one_(mat_T) == FALSE)) {
    stop("hmm_simulate_normal_data(): rows in the transition matrix `mat_T` must sum up to 1", call. = FALSE)
  }
  
  if (length(L) != 1) {
    stop("hmm_simulate_normal_data(): the number of simulations `L` must be a single integer", call. = FALSE)
  }
  
  vstate <- rep(NA, L)
  vemit <- rep(NA, L)
  
  n_states <- nrow(mat_T)
  pi <- get_pi_(mat_T)
  
  vstate[1] <- sample.int(n = n_states, size = 1, replace = TRUE, prob = pi)
  vemit[1] <- stats::rnorm(n = 1, mean = means[vstate[1]], sd = sigma)
  
  for(i in 2:L) {
    p_state <- mat_T[vstate[i - 1], ]
    vstate[i] <- sample.int(n = n_states, size = 1, replace = TRUE, prob = p_state)
    vemit[i] <- stats::rnorm(n = 1, mean = means[vstate[i]], sd = sigma)
  }
  
  list("data" = vemit,
       "states" = factor(vstate, levels = 1:n_states),
       "pi" = pi)
}


#' Forward-backward algorithm to calculate the posterior probabilities of hidden states.
#'
#' Forward-backward algorithm to calculate the posterior probabilities of hidden states.
#'
#' @param data (numeric) normal data
#'
#' @param pi (numeric) prior probability of states
#'
#' @param mat_T (matrix) transition probability matrix
#'
#' @param means (numeric) vector with prior means
#'
#' @param sdev (numeric) prior standard deviation
#'
#' @details
#' Here details on how the calculation is made
#'
#' @references
#' Here some references
#'
#' @return
#' List with posterior probabilities (TO BE CORRECTED)
#'
#' @export
#'
#' @examples
#' prior_mat <- rbind(c(1-0.05, 0.05, 0),
#'                   c(0.05, 1-0.1, 0.05),
#'                   c(0, 0.05, 1-0.05))
#'
#' prior_means <- c(-0.1, 0.0, 0.1)
#' prior_sd  <- sqrt(0.1)
#' L <- 100
#'
#' # Simulate HMM model based on normal data based on prior information
#' sim_data_normal <- hmm_simulate_normal_data(L = L,
#'                                             mat_T = prior_mat,
#'                                             means = prior_means,
#'                                             sigma = prior_sd)
#' pi <- sim_data_normal$pi
#' # pi <- get_pi(prior_mat)
#' hmm_norm_data <- sim_data_normal$data
#'
#' # Calculate posterior probabilities of hidden states
#' post_prob <-  posterior_probabilities_normal(data = hmm_norm_data,
#'                                              pi = pi,
#'                                              mat_T = prior_mat,
#'                                              means = prior_means,
#'                                              sdev = prior_sd)

posterior_probabilities_normal <- function(data, pi, mat_T, means, sdev) {
  
  cap_ <- .Machine$double.xmax
  floor_ <- .Machine$double.xmin
  
  L <- length(data)
  n_states <- nrow(mat_T)
  mat_F <- matrix(ncol = n_states, nrow = L)
  vs <- rep(1, L) # vector for scaling
  
  pi[pi < floor_] <- floor_ # never 0 so that there is no 0*Inf. 1e-100*Inf = Inf
  vtemp <- cap_floor_(pi * stats::dnorm(x = data[1], mean = means, sd = sdev), cap_, floor_)
  vs[1] <- sum(vtemp) # scaling: Pr(y[1]), Pr(y[2]|y[1]),..., Pr(y[i+1]|y[i]), ...
  mat_F[1,] <- vtemp / vs[1]
  
  # Forward
  for (l in 2:L) {
    probs <- cap_floor_(stats::dnorm(x = data[l], mean = means, sd = sdev), cap_, floor_)
    vtemp <- mat_F[l - 1,] %*% mat_T * probs
    vs[l] <- sum(vtemp)
    mat_F[l,] <- vtemp / vs[l]
  }
  
  # Backward
  mat_B <- matrix(nrow = L, ncol = n_states)
  mat_B[L, ] <- rep(1, n_states) / vs[L]
  
  for (l in L:2) {
    probs <- cap_floor_(stats::dnorm(x = data[l], mean = means, sd = sdev), cap_, floor_)
    mat_B[l - 1, ] <- mat_T %*% (probs * mat_B[l, ]) / vs[l - 1]
  }
  
  list(F = mat_F, B = mat_B, s = vs)
}


# sample_states_normal_ - changes:
# - the function returns states as factor (helps to calculate means and variances for missing states as well as calculating table with state transitions)
# - is made more robust (no zero probabilities are produced)

##!!! LCM: Should warnings be given in the function below when probabilities are padded?

# MM:
# 1) states are returned as factors (factor(states, levels = 1:DIM))
#    table(states) always includes all states, even if they are not found in
#    the vector (levels are stored as meta-attributes). Example below:
#      x <- factor(1:10)
#      attributes(x)
#      table(x[1:2])
# 2) dnorm is used in a vectorized way - probabilities are calculated only once
# 3) sample() is replaced with two-three fold more efficient alternative
#
# 4) ad.3) the alternative is further improved in terms of efficiency:
#     - instead of generating L uniforms, all uniforms are generated at once,
#       stored in a variable, which is then subsetted in each for-loop iteration.
#       --> almost two-fold speed improvement (~1.7x faster for L=1024 and L=32768)
#
#    all_runif <- runif(L)
#    states[1] <- which.max(cumsum(p) - all_runif[1] > 0)
#    for (...)
#    states[l] <- which.max(cumsum(p) - stats::runif(1) > 0)

#' @keywords internal
sample_states_normal_ <- function(mat_R, mat_T, pi, means, sdev, data) {
  
  L <- length(mat_R$F[ ,1])
  n_states <- length(mat_R$F[1, ])
  states <- rep(NA, L)
  cap_ <- .Machine$double.xmax
  floor_ <- .Machine$double.xmin
  
  p <- cap_floor_(mat_R$F[1, ] * mat_R$F[1, ] * mat_R$s[1], cap_, floor_)
  p <- p / sum(p)
  p[is.na(p)] <- floor_
  all_runif <- stats::runif(L) # generating L unifs at once is faster than generating L times 1 unif
  states[1] <- which.max(cumsum(p) - all_runif[1] > 0)
  # states[1] <- which.max(cumsum(p) - stats::runif(1) > 0)
  
  mat_means <- matrix(means, nrow = L, ncol = length(means), byrow = TRUE)
  probs <- cap_floor_(stats::dnorm(data, mat_means, sdev), cap_, floor_)
  
  for (l in 2:L) {
    p <- mat_T[states[l - 1], ] * probs[l, ] * mat_R$B[l, ]
    p <- p / sum(p)
    p[is.na(p)] <- floor_
    states[l] <- which.max(cumsum(p) - all_runif[l] > 0)
    # states[l] <- which.max(cumsum(p) - stats::runif(1) > 0)
  }
  factor(states, levels = 1:n_states)
}

# sample_T_ - changes:
# - shorter code but not necessarily faster
# CV: function should return a sample from the prior if states==NULL
# !!! LCM: Should warnings be given in the function below when probabilities are padded?

#' @keywords internal
sample_T_ <- function(prior_mat, states = NULL) {
  
  n_states <- length(prior_mat[1, ])
  mat_counts <- prior_mat
  
  if (!is.null(states)) {
    L <- length(states)
    mat_counts <- mat_counts + table(states[1:(L - 1)], states[2:L])
  }
  
  # Extract the diagonal
  vtemp <- diag(mat_counts, names = FALSE)
  
  # Extract the off-diagonal
  for (i in 2:n_states) {
    vtemp[n_states + i - 1] <- mat_counts[i-1,i] + mat_counts[i,i - 1]
  }
  
  floor_ <- .Machine$double.xmin
  vtemp[vtemp == 0] <- floor_
  vtemp[is.na(vtemp)] <- floor_
  
  # vtemp <- MCMCprecision::rdirichlet(n = 1, a = vtemp)
  vtemp <- gtools::rdirichlet(n = 1, alpha = vtemp)
  mat_T <- diag(vtemp[1:n_states])
  
  for (i in 2:n_states) {
    mat_T[i,i-1] <- vtemp[n_states + i - 1] * 0.5
    mat_T[i-1,i] <- vtemp[n_states + i - 1] * 0.5
  }
  
  for (i in 1:n_states) {
    mat_T[i, ] <- mat_T[i, ] / sum(mat_T[i, ])
  }
  mat_T
}

# sample_means_sd_ - changes:
# - recovering sums of squares out of variance with correct multiplicative factor (n-1)
# - when some state is not represented (count=0) then the prior parameters are used only
# - function is 2x-3x faster than the original one
# CV: this function needs to work also without data and then return a sample from the prior

#' @keywords internal

# LCM: BELOW IS THE NEW SAMPLING FUNCTION FOR VERSION 8
sample_means_sd_ <- function(prior_means, prior_var, states = NULL, data = NULL) {
  L <- length(states)
  n_states <- length(prior_means)
  ss <- 0 
  ybar <- 0
  n <- 0
  if (!is.null(states) & !is.null(data)) {
    n <- tapply(data, states, length, default = 0)
    ybar <- tapply(data, states, mean, default = 0)
    vars <- tapply(data, states, stats::var, default = 0) # Denominator n-1
    vars[n == 1] <- 0 
    ss <- sum(vars * (n - 1)) ## NEW V9: since the denominator is n-1 should not the factor here be n-1????
  }      
  # Inverse gamma with shape=(L+1)/2, scale=(prior_var+ss)/2   
  # CV: I cannot see where the variable "L" is defined; (FIXED) 
  #note: sigma2 is the variance!
  sigma2 <- 1 / stats::rgamma(n = 1, shape = (L + 1) * 0.5,
                              rate = (prior_var + ss + sum((ybar-prior_means)*(ybar-prior_means)*n/(n+1))) * 0.5) 
  means <- stats::rnorm(n = n_states, 
                        mean = (n * ybar  + prior_means) / (n + 1),
                        ## sd = sigma2 / (L + 1)) ## BUGBUG: this is the variance and not the sdev!
                        sd = sqrt(sigma2/(n+1))) 
  list(means = sort(means), sdev = sqrt(sigma2))
}




#' @keywords internal
#### NEW LCM: remove below????
#get_parms_log_prob_ <- function(means, sd, mat_T, prior_means, prior_var, prior_T) {

#    lp <- stats::dgamma(1 / sd^2, shape = 1 * 0.5, rate = prior_var * 0.5)
#    lp <- lp + sum(stats::dnorm(x = (means - prior_means) / sd, sd = sd, log = TRUE))
#    lp <- lp + gtools::ddirichlet(x = mat_T[1,1:2], alpha = prior_T[1,1:2])
#    K <- length(prior_means)

#    if (K > 2) {
#        for (i in 2:(K-1)) {
#            lp <- lp + gtools::ddirichlet(x = mat_T[i,(i-1):(i+1)], alpha = prior_T[i,(i-1):(i+1)])
#        }
#    }
#    lp <- lp + gtools::ddirichlet(x = mat_T[K,(K-1):K], alpha = prior_T[K,(K-1):K])
#    lp
#}


#' @keywords internal
init_hmm_mcmc_normal_ <- function(data, prior_T, prior_means, prior_sd,
                                  init_T=NULL, init_means=NULL, init_sd=NULL, verbose,
                                  iter, warmup, thin, chain_id = NULL) {
  ## CV: init_T, init_means and init_sd may be set to NULL, 
  ## in this case the initial values are sampled from the prior 
  
  if (!is.numeric(data)) {
    stop("hmm_mcmc_normal(): `data` needs to be a numeric vector", call. = FALSE)
  }
  ## new LCM
  if (sum(is.na(data))>0){
    stop("hmm_mcmc_normal(): `data` contains missing values", call. = FALSE)
  }
  ##end  new LCM
  if (any(is_row_sum_one_(prior_T) == FALSE)) {
    stop("hmm_mcmc_normal(): rows in the transition matrix `prior_T` must sum up to 1", call. = FALSE)
  }
  if (length(prior_means) != nrow(prior_T) | length(prior_means) != ncol(prior_T)) {
    stop("hmm_mcmc_normal(): number of states is not the same between input variables", call. = FALSE)
  }
  if (length(prior_sd) != 1) {
    stop("hmm_mcmc_normal(): initial standard deviation `prior_sd` must be of length 1", call. = FALSE)
  }
  if (warmup >= iter + 2) {
    stop("hmm_mcmc_normal(): `warmup` must be lower than `iter`", call. = FALSE)
  }
  if (thin > iter - warmup) {
    stop("hmm_mcmc_normal(): `thin` cannot exceed iterations after warmup period", call. = FALSE)
  }
  
  # Initialization
  # CV: function returns a sample from the prior
  out <- sample_means_sd_(prior_means = prior_means, prior_var = prior_sd^2)
  
  if (is.null(init_means)) {
    init_means <- out$means
  }
  if (is.null(init_sd)) {
    init_sd <- out$sdev
  }
  if (is.null(init_T)) {
    init_T <- sample_T_(prior_mat = prior_T) ## CV: function returns a sample from the prior
  }
  if (length(init_means) != nrow(init_T) | length(init_means) != ncol(init_T)) {
    stop("hmm_mcmc_normal(): number of states is not the same between input variables", call. = FALSE)
  }
  if (length(init_sd) != 1) {
    stop("hmm_mcmc_normal(): initial standard deviation `init_sd` must be of length 1", call. = FALSE)
  }
  if (any(is_row_sum_one_(init_T) == FALSE)) {
    stop("hmm_mcmc_normal(): rows in the transition matrix `init_T` must sum up to 1", call. = FALSE)
  }
  
  # TODO: estimate more or less how far can prior mean be
  # from the sample mean such that the density is still > 0
  sample_mean <- mean(data)
  sample_sd <- stats::sd(data)
  abs_mean_ratios <- abs(init_means) / abs(sample_mean)
  sd_ratio <- init_sd / sample_sd
  
  chain_char <- NULL
  if (!is.null(chain_id)) {
    chain_char <- paste0("(chain ", chain_id, ") ")
  }
  if (verbose) {
    # MM:
    # CHANGED: warnings --> message
    # Reason: message is outputted to the console immediately
    #         warnings are given after the function finishes all calculations
    
    if (any(abs_mean_ratios > 5)) {
      message("hmm_mcmc_normal(): ", chain_char, "at least one element in `init_means` is at least 5x bigger than the sample mean.")
    }
    if (sd_ratio > 5) {
      message("hmm_mcmc_normal(): ", chain_char, "`init_sd` at least 5x bigger than the sample standard deviation.")
    }
  }
  
  init_pi <- get_pi_(init_T)
  
  init_mat_res <- posterior_probabilities_normal(data = data,
                                                 pi = init_pi,
                                                 mat_T = init_T,
                                                 means = init_means,
                                                 sdev = init_sd)
  init_states <- sample_states_normal_(mat_R = init_mat_res,
                                       pi = init_pi,
                                       mat_T = init_T,
                                       means = init_means,
                                       sdev = init_sd,
                                       data = data)
  
  list(init_states = init_states, init_means = init_means, init_sd = init_sd,
       init_T = init_T, init_pi = init_pi, init_mat_res = init_mat_res)
}


#' MCMC simulation of a Hidden Markov Normal Model
#'
#'
#' @param data (numeric) normal data
#'
#' @param prior_T (matrix) prior transition matrix
#'
#' @param prior_means (numeric) prior means
#'
#' @param prior_sd (numeric) a single prior standard devation
#'
#' @param iter (integer) number of MCMC iterations
#'
#' @param warmup (integer) number of warmup iterations
#'
#' @param thin (integer) thinning parameter. By default, \code{1}
#'
#' @param seed (integer) seed parameter
#'
#' @param init_T (matrix) \code{optional parameter}; initial transition matrix
#'
#' @param init_means (numeric) \code{optional parameter}; initial means
#'
#' @param init_sd (numeric) \code{optional parameter}; initial standard deviation
#'
#' @param print_params (logical) \code{optional parameter}; print parameters every iteration. By default, \code{TRUE}
#'
#' @param verbose (logical) \code{optional parameter}; print additional messages. By default, \code{TRUE}
#'
#' @details
#' Here details
#'
#' @references
#' Here references
#'
#' @return
#' List with following elements:
#' \itemize{
#'   \item data: data used for simulation
#'   \item estimates: list with various estimates
#'   \item idx: indices with iterations after the warmup period
#'   \item priors: prior parameters
#'   \item inits: initial parameters
#'   \item last_iter: list with samples from the last MCMC iteration
#'   \item info: list with various meta information about the object
#' }
#'
#' @export
#'
#' @examples
#' # Simulate normal data
#' N <- 2^10
#' true_T <- rbind(c(0.95, 0.05, 0),
#'                 c(0.025, 0.95, 0.025),
#'                 c(0.0, 0.05, 0.95))
#'
#' true_means <- c(-5, 0, 5)
#' true_sd <- 1.5
#'
#' simdata_full <- hmm_simulate_normal_data(L = N, 
#'                                          mat_T = true_T, 
#'                                          means = true_means,
#'                                          sigma = true_sd)
#' simdata <- simdata_full$data
#' hist(simdata, breaks = 40, probability = TRUE,  
#'      main = "Distribution of the simulated normal data")
#' lines(density(simdata), col = "red")
#'
#' # Set numbers of states to be inferred
#' n_states_inferred <- 3
#' 
#' # Set priors
#' prior_T <- generate_random_T(n_states_inferred)
#' prior_means <- c(-18, -1, 12)
#' prior_sd <- 3
#' 
#' # Simmulation settings
#' iter <- 50
#' warmup <- floor(iter / 5) # 20 percent
#' thin <- 1
#' seed <- sample.int(10000, 1)
#' print_params <- FALSE # if TRUE then parameters are printed in each iteration
#' verbose <- FALSE # if TRUE then the state of the simulation is printed
#' 
#' # Run MCMC sampler
#' res <- hmm_mcmc_normal(data = simdata,
#'                        prior_T = prior_T,
#'                        prior_means = prior_means,
#'                        prior_sd = prior_sd,
#'                        iter = iter,
#'                        warmup = warmup,
#'                        seed = seed,
#'                        print_params = print_params,
#'                        verbose = verbose)
#' res
#' 
#' summary(res) # summary output can be also assigned to a variable
#' 
#' coef(res) # extract model estimates
#' 
#' # plot(res) # MCMC diagnostics

hmm_mcmc_normal <- function(data,
                            prior_T,
                            prior_means,
                            prior_sd,
                            iter = 600,
                            warmup = floor(iter / 5),
                            thin = 1,
                            seed = sample.int(.Machine$integer.max, 1),
                            init_T = NULL,
                            init_means = NULL,
                            init_sd = NULL,
                            print_params = TRUE,
                            verbose = TRUE) {
  
  set.seed(seed)
  
  init_data <- init_hmm_mcmc_normal_(data, prior_T, prior_means, prior_sd,
                                     init_T, init_means, init_sd, verbose,
                                     iter, warmup, thin)
  
  n_data <- length(data)
  n_states <- length(prior_means)
  
  mean_states <- matrix(0, ncol = n_states, nrow = n_data)
  all_means <- matrix(0, ncol = n_states, nrow = iter)
  all_sd <- numeric(iter)
  all_mat_T <- array(dim = c(dim(prior_T), iter))
  vlh <- numeric(iter)
  
  all_means[1, ] <- init_data$init_means
  all_sd[1] <- init_data$init_sd
  all_mat_T[, ,1] <- init_data$init_T
  vlh[1] <- sum(log(init_data$init_mat_res$s))
  states <- init_data$init_states
  prior_pi=get_pi_(prior_T) ## New V9
  prior_P=diag(prior_pi)%*%prior_T ## New V9
  
  # Run sampler
  for (it in 2:iter) {
    
    print_progress_(it, iter, verbose = verbose)
    
    m <- sample_means_sd_(states = states,
                          prior_means = prior_means,
                          prior_var = prior_sd^2,
                          data = data)
    means <- m$means
    sd <- m$sdev
    mat_T <- sample_T_(states, prior_mat = prior_P)
    pi <- get_pi_(mat_T)
    mat_res <- posterior_probabilities_normal(data = data,
                                              pi = pi,
                                              mat_T = mat_T,
                                              means = means,
                                              sdev = sd)
    states <- sample_states_normal_(mat_R = mat_res,
                                    pi = pi,
                                    mat_T = mat_T,
                                    means = means,
                                    sdev = sd,
                                    data = data)
    if (iter > warmup) {
      for (i in 1:n_states) {
        mean_states[ ,i] <- mean_states[ ,i] + (states == i)
      }
    }
    
    all_means[it, ] <- means
    all_sd[it] <- sd
    all_mat_T[ , ,it] <- mat_T
    vlh[it] <- sum(log(mat_res$s))
    
    if (print_params) {
      print(c(colMeans(all_means[1:it, ]), mean(all_sd[it])))
    }
  }
  
  # Prepare outputs
  idx <- seq.int(warmup + 1, by = thin, to = iter)
  
  #### NEW LCM: remove below????
  #vlparms <- numeric(length(idx))
  #for(i in 1:length(idx)) {
  #    vlparms[i] <- get_parms_log_prob_(means = all_means[warmup + i],
  #                                      sd = all_sd[warmup + i],
  #                                      mat_T = all_mat_T[, , warmup + i],
  #                                      prior_means = prior_means,
  #                                      prior_var = (prior_sd)^2,
  #                                      prior_T = prior_T)
  #}
  
  colnames(all_means) <- paste0("mean[", 1:n_states,"]")
  samples <- list(means = all_means,
                  sd = all_sd,
                  mat_T = all_mat_T,
                  vlh = vlh)
  
  posterior_states <- factor(max.col(mean_states, "first"), 1:n_states)
  
  estimates <- list(means = colMeans(all_means[idx, ]),
                    sd = mean(all_sd[idx]),
                    mat_T = apply(all_mat_T[ , ,idx], c(1,2), mean),
                    posterior_states = posterior_states,
                    posterior_states_prob = mean_states / iter,
                    #log_posterior = vlparms, # NEW LCM: REMOVED AGAIN _FOR BAYES INTERPRETATION
                    log_likelihood = vlh[idx])    
  
  priors <- list(prior_means = prior_means,
                 prior_sd = prior_sd,
                 prior_T = prior_T)
  
  inits <- list(init_states = init_data$init_states,
                init_means = init_data$init_means,
                init_sd = init_data$init_sd,
                init_T = init_data$init_T)
  
  last_iter <- list(means = all_means[iter, ],
                    sd = all_sd[iter],
                    mat_T = all_mat_T[ , ,iter])
  
  info <- list(model_name = "hmm_mcmc_normal",
               date = as.character(Sys.time()),
               seed = seed,
               iter = iter,
               warmup = warmup,
               thin = thin)
  
  res <- list(data = data,
              samples = samples,
              estimates = estimates,
              idx = idx,
              priors = priors,
              inits = inits,
              last_iter = last_iter,
              info = info)
  
  class(res) <- "hmm_mcmc_normal"
  res$info$object_size <- format(utils::object.size(res), "Mb")
  res
}


#' @keywords internal
#' @export
print.hmm_mcmc_normal <- function(x, ...) {
  info <- x$info
  mod <- if (info$model_name == "hmm_mcmc_normal") "HMM Normal" else NA
  cat("Model:", mod, "\n")
  cat("Type:", "MCMC", "\n")
  cat("Iter:", info$iter, "\n")
  cat("Warmup:", info$warmup, "\n")
  cat("Thin:", info$thin, "\n")
  cat("States:", length(x$priors$prior_means), "\n")
}


#' @keywords internal
#' @export
summary.hmm_mcmc_normal <- function(object, ...) {
  
  info <- object$info
  idx <- object$idx
  data <- object$data
  
  m_est <- object$estimates$means
  sd_est <- object$estimates$sd
  T_est <- object$estimates$mat_T
  
  post_states <- object$estimates$posterior_states
  state_tab <- table(post_states, dnn = "")
  
  dens_data <- stats::density(data)
  kl_list <- rep(NA, 500)
  for (j in 1:500) {
    sim_output <- unlist(lapply(1:length(state_tab), function(i) {
      stats::rnorm(state_tab[i], m_est[i], sd_est) })
    )
    dens_sim <- stats::density(sim_output)
    kl_list[j] <- kullback_leibler_cont_appr(dens_data$y, dens_sim$y)
  }
  kl_div <- mean(kl_list)
  
  ll_info <- c(mean(object$estimates$log_likelihood),
               stats::sd(object$estimates$log_likelihood),
               stats::median(object$estimates$log_likelihood))
  names(ll_info) <- c("mean", "sd", "median")
  
  #### NEW LCM: remove below????
  #post_info <- c(mean(object$estimates$log_posterior),
  #               stats::sd(object$estimates$log_posterior),
  #               stats::median(object$estimates$log_posterior))
  # names(post_info) <- c("mean", "sd", "median")
  
  summary_res <- list("estimated_means" = m_est,
                      "estimated_sd" = sd_est,
                      "estimated_transition_rates" = T_est,
                      "assigned_states" = state_tab,
                      "approximate_kullback_leibler_divergence" = kl_div,
                      "log_likelihood" = ll_info)
  #### NEW LCM: remove below????
  #"log_posterior" = post_info)
  
  cat("Estimated means:\n")
  print(summary_res$estimated_means)
  cat("\n")
  
  cat("Estimated sd:\n")
  cat(summary_res$estimated_sd)
  cat("\n")
  cat("\n")
  
  cat("Estimated transition rates:\n")
  etr <- summary_res$estimated_transition_rates
  rownames(etr) <- colnames(etr) <- names(summary_res$assigned_states)
  print(etr)
  cat("\n")
  
  cat("Assigned states:\n")
  as <- summary_res$assigned_states
  as_names <- attributes(as)$dimnames[[1]]
  print(stats::setNames(as.numeric(as), as_names))
  cat("\n")
  
  cat("Approximate Kullback-Leibler divergence:\n")
  cat(stats::setNames(summary_res$approximate_kullback_leibler_divergence, ""))
  cat("\n")
  cat("\n")
  
  cat("Log Likelihood:\n")
  print(summary_res$log_likelihood)
  cat("\n")
  
  #### NEW LCM: remove below????
  #cat("Log Posterior:\n")
  #print(summary_res$log_posterior)
  
  invisible(summary_res)
}


#' Extract model estimates
#'
#' \code{coef} is a generic function which extracts model estimates from \code{mcmc_hmm_*} objects
#'
#' @param object an object of class inheriting from "\code{mcmc_hmm_*}"
#'
#' @param ... not used
#'
#' @return Estimates extracted from MCMC HMM objects
#'
#' @export
#' @export coef.hmm_mcmc_normal
#'
#' @examples
#' # TODO

coef.hmm_mcmc_normal <- function(object, ...) {
  est <- object$estimates
  list(means = est$means,
       sd = est$sd,
       mat_T = est$mat_T)
}


### NOTE: Can one suppress the counts this function spits out while running/plotting the list of plots?


#' Plot method for \code{hmm_mcmc_normal} objects
#'
#' @param x (hmm_mcmc_\*) MCMC HMM object
#'
#' @param simulation (logical)
#'
#' @param true_means (numeric)
#'
#' @param true_sd (numeric)
#'
#' @param true_mat_T (matrix)
#'
#' @param true_states (integer)
#'
#' @param ... not used
#'
#' @details
#' Here details
#'
#' @return
#' No return value
#'
#' @export
#' @export plot.hmm_mcmc_normal
#'
#' @importFrom stats qqplot
#'
#' @examples
#' # TODO
#'

# simulation = TRUE
# true_T <- rbind(c(0.95, 0.05, 0),
#                 c(0.025, 0.95, 0.025),
#                 c(0.0, 0.05, 0.95))
#
# true_means <- c(-5, 0, 5)
# true_sd <- 1.5
# true_states <- simdata_full$states

plot.hmm_mcmc_normal <- function(x,
                                 simulation = FALSE,
                                 true_means = NULL,
                                 true_sd = NULL,
                                 true_mat_T = NULL,
                                 true_states = NULL,
                                 ...) {
  
  info <- x$info
  data <- x$data
  idx <- x$idx
  
  if (simulation) {
    cond <- any(c(is.null(true_means), is.null(true_sd), is.null(true_mat_T),
                  is.null(true_states)))
    if (cond) {
      stop("plot.hmm_mcmc_normal(): if `simulation=TRUE` then `true_means` ",
           "`true_sd`, `true_mat_T` and `true_states` must be defined", call. = FALSE)
    }
  }
  
  # Diagnostics mean
  all_means <- convert_to_ggmcmc(x, pattern = "mean")
  n_means <- attributes(all_means)$nParameters
  facet_means <- ggplot2::facet_wrap(~ Parameter, ncol = floor(n_means / 2), scales = "free")
  mdens <- ggmcmc::ggs_density(all_means) + facet_means
  mtrace <- ggmcmc::ggs_traceplot(all_means) + facet_means
  mauto <- invisible(utils::capture.output(ggmcmc::ggs_autocorrelation(all_means))) +
  facet_means
  
  # Diagnostics transitions
  all_T <- convert_to_ggmcmc(x, pattern = "T")
  n_t <- attributes(all_T)$nParameter
  facet_t <- ggplot2::facet_wrap(~ Parameter, ncol = sqrt(n_t), scales = "free")
  labels_t <- ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = "."))
  Ttrace <- ggmcmc::ggs_traceplot(all_T) + facet_t + labels_t + ggplot2::labs(x = "Iteration", y = "Value")
  Tdens <- ggmcmc::ggs_density(all_T) + facet_t + labels_t + ggplot2::labs(x = "Value", y = "Density")
  
  # Diagnostics sd
  df_sigma <- convert_to_ggmcmc(x, "sigma")
  sdtrace <- ggmcmc::ggs_traceplot(df_sigma) + ggplot2::labs(x = "Iteration", y = "Value")
  sddens <- ggmcmc::ggs_density(df_sigma) + ggplot2::labs(x = "Value", y = "Density")
  
  # Likelihood trace
  lltrace <- x$estimates$log_likelihood
  lltrace_df <- as.data.frame(cbind(c((info$warmup+1):info$iter),lltrace))
  names(lltrace_df) <- c("iteration", "log_likelihood")
  
  # llplot <- ggplot2::ggplot(lltrace_df, ggplot2::aes(x = iteration, y = log_likelihood)) +
  #   ggplot2::geom_line() +
  #   ggplot2::labs(x = "Iteration", y = "Log-likelihood")
  
  llplot <- ggplot2::ggplot(lltrace_df, ggplot2::aes_string(x = "iteration", y = "log_likelihood")) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iteration", y = "Log-likelihood")
  
  # NEW LCM: REMOVE BELOW
  # lposterior trace
  #posttrace <- x$estimates$log_posterior
  #posttrace_df <- as.data.frame(cbind(idx,posttrace))
  #names(posttrace_df) <- c("iteration", "log_posterior")
  
  #postplot <- ggplot2::ggplot(posttrace_df, ggplot2::aes(x = iteration, y = log_posterior)) +
  #    ggplot2::geom_line() +
  #    ggplot2::labs(x = "Iteration", y = "Log-posterior")
  
  # Confusion matrix
  if (simulation) {
    m_multi <- list("target" = true_states,
                    "prediction" = x$estimates$posterior_states)
    conf_mat <- cvms::confusion_matrix(targets = m_multi$target,
                                       predictions = m_multi$prediction)
    
    conf_mat_plot <- suppressWarnings(
      cvms::plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]],
                                  add_sums = TRUE)
    )
  }
  
  # Check assignment of states along chromosome
  states_df <- as.data.frame(cbind(1:length(data), data, x$estimates$posterior_states))
  post_means <- numeric(length(data))
  
  for (l in 1:length(data)) {
    post_means[l] <- sum(x$estimates$means * x$estimates$posterior_states_prob[l, ])
  }
  
  states_df$post_means <- post_means
  names(states_df) <- c("position", "data", "posterior_states", "posterior_means")
  states_df$posterior_states <- as.factor(states_df$posterior_states)
  
  # statesplot <- ggplot2::ggplot(states_df, ggplot2::aes(x = position, y = data)) +
  #   ggplot2::geom_line(col = "grey") +
  #   ggplot2::geom_point(ggplot2::aes(colour = posterior_states), shape = 20, size = 1.5, alpha = 0.75) +
  #   ggplot2::geom_line(ggplot2::aes(x = position, y = posterior_means), size = 0.15) +
  #   ggplot2::guides(colour = ggplot2::guide_legend(title = "Post States")) +
  #   ggplot2::labs(x = "Position", y = "Data")
  
  statesplot <- ggplot2::ggplot(states_df, ggplot2::aes_string(x = "position", y = "data")) +
    ggplot2::geom_line(col = "grey") +
    ggplot2::geom_point(ggplot2::aes_string(colour = "posterior_states"), shape = 20, size = 1.5, alpha = 0.75) +
    ggplot2::geom_line(ggplot2::aes_string(x = "position", y = "posterior_means"), size = 0.15) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Post States")) +
    ggplot2::labs(x = "Position", y = "Data")
  
  if (simulation) {
    states_df2 <- as.data.frame(cbind(1:length(data), data, true_states))
    post_means <- true_means[true_states]
    states_df2$post_means <- post_means
    names(states_df2) <- c("position", "data", "true_states", "posterior_means")
    states_df2$true_states <- as.factor(states_df2$true_states)
    
    # statesplot2 <- ggplot2::ggplot(states_df2, ggplot2::aes(x = position, y = data)) +
    #   ggplot2::geom_line(col = "grey") +
    #   ggplot2::geom_point(ggplot2::aes(colour = true_states), shape = 20, size = 1.5, alpha = 0.75) +
    #   ggplot2::geom_line(ggplot2::aes(x = position, y = posterior_means), size = 0.15) +
    #   ggplot2::guides(colour = ggplot2::guide_legend(title = "True States")) +
    #   ggplot2::labs(x = "Position", y = "Data")
    
    statesplot2 <- ggplot2::ggplot(states_df2, ggplot2::aes_string(x = "position", y = "data")) +
      ggplot2::geom_line(col = "grey") +
      ggplot2::geom_point(ggplot2::aes_string(colour = "true_states"), shape = 20, size = 1.5, alpha = 0.75) +
      ggplot2::geom_line(ggplot2::aes_string(x = "position", y = "posterior_means"), size = 0.15) +
      ggplot2::guides(colour = ggplot2::guide_legend(title = "True States")) +
      ggplot2::labs(x = "Position", y = "Data")
  }
  
  # qq-plot
  dims <- length(x$estimates$means)
  dist <- rep("norm", dims)    #### NOTE: distribution here will have to be changed for count version!!
  params <- list()
  for (j in 1:dims) {
    params[[j]] <- c(as.numeric(x$estimates$means[j]), x$estimates$sd)
  }
  tab_post <- table(x$estimates$posterior_states)
  weight <- as.numeric(tab_post / sum(tab_post))
  resulting_mixture_of_normals <- mistr::mixdist(dist, params, weights = weight)
  
  qqplot <- mistr::QQplotgg(data, resulting_mixture_of_normals, col = "black", line_col = "blue") +
    ggplot2::labs(x = "Resulting mixture of normals", y = "Data") +
    ggplot2::theme_get()
  
  # Compare densities
  # Note: both rnorm() and density() must be adapted in count version!!!
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
  dens_df$value <- as.numeric(dens_df$value)
  # kl_plot <- ggplot2::ggplot(dens_df, ggplot2::aes(x = as.numeric(value), fill = data_type)) +
  #   ggplot2::geom_density(alpha = 0.4) +
  #   ggplot2::labs(title = "Model Fit", x = "Values", y = "Density") +
  #   ggplot2::geom_vline(xintercept = x$estimates$means, color = "red", size = 1) +
  #   ggplot2::geom_vline(xintercept = c(x$estimates$means) + x$estimates$sd,
  #                       linetype = "dotted", color = "red", size = 1) +
  #   ggplot2::geom_vline(xintercept = c(x$estimates$means) - x$estimates$sd,
  #                       linetype = "dotted",  color = "red", size = 1) +
  #   ggplot2::guides(fill = ggplot2::guide_legend(title = "Data type"))
  
  
  kl_plot <- ggplot2::ggplot(dens_df, ggplot2::aes_string(x = "value", fill = "data_type")) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::labs(title = "Model Fit", x = "Values", y = "Density") +
    ggplot2::geom_vline(xintercept = x$estimates$means, color = "red", size = 1) +
    ggplot2::geom_vline(xintercept = c(x$estimates$means) + x$estimates$sd,
                        linetype = "dotted", color = "red", size = 1) +
    ggplot2::geom_vline(xintercept = c(x$estimates$means) - x$estimates$sd,
                        linetype = "dotted",  color = "red", size = 1) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Data type"))
  
  if (simulation) {
    kl_plot <- kl_plot +
      ggplot2::geom_vline(xintercept = true_means, color = "blue", size = 0.5) +
      ggplot2::geom_vline(xintercept = c(true_means) + true_sd,  linetype = "dotted",
                          color = "blue", size = 0.5) +
      ggplot2::geom_vline(xintercept = c(true_means) - true_sd, linetype = "dotted",
                          color = "blue", size = 0.5)
  }
  
  if (simulation) {
    plotlist <- list(mtrace, mdens, mauto, Ttrace, Tdens, sdtrace, sddens,
                     llplot, conf_mat_plot, qqplot, kl_plot)
  } else {
    plotlist <- list(mtrace, mdens, mauto, Ttrace, Tdens, sdtrace, sddens,
                     llplot, statesplot, qqplot, kl_plot)
  }
  
  
  
  for (ii in 1:length(plotlist)) {
    
    if (simulation && ii == 10) {
      gridExtra::grid.arrange(statesplot, statesplot2, nrow = 2)
    }
    invisible(utils::capture.output(print(plotlist[[ii]])))
  }
}
# plot(res)
# plot(res, simulation = TRUE, true_means = true_means,
#      true_sd = true_sd, true_mat_T = true_mat_T,
#      true_states = true_states)


#' Calculate a confusion matrix...DESCRIPTION TO BE IMPROVED
#'
#' A diagnostic function that tests the reliability of estimation
#' procedures given the inferred transition rates
#'
#' @param N (numeric) number of simulations
#'
#' @param res (mcmc_hmm_\*) simulated MCMC HMM model
#'
#' @param plot (logical) plot confusion matrix. By default \code{TRUE}
#'
#' @details
#' First the data is simulated given the inferred model parameters and transition
#' rates. Then posterior probabilities are calculated and states are inferred.
#' Finally, the inferred states and simulated states are compared via
#' \code{\link[cvms]{confusion_matrix}} function.
#'
#' @return
#' Confusion matrix: \code{\link[cvms]{confusion_matrix}}
#'
#' @export
#'
#'
#' @examples
#' # TODO

conf_mat <- function(N, res, plot = TRUE) {
  
  if (inherits(res, "hmm_mcmc_normal")) {
    trial <- hmm_simulate_normal_data(N,
                                      res$estimates$mat_T,
                                      res$estimates$means,
                                      res$estimates$sd)
  } else if (inherits(res, "hmm_mcmc_poisson")) {
    trial <- hmm_simulate_poisgamma_data(N,
                                         res$estimates$mat_T,
                                         res$estimates$betas,
                                         res$estimates$alpha)
  } else {
    stop("conf_mat(): currently \"hmm_mcmc_normal\" and \"hmm_mcmc_poisson\" models are supported", call = FALSE)
  }
  
  vpi <- trial$pi
  trial_data <- trial$data
  trial_states <- trial$states
  
  if (inherits(res, "hmm_mcmc_normal")) {
    # Calculate posterior probabilities of hidden states
    post_prob <-  posterior_probabilities_normal(data = trial_data,
                                                 pi = vpi,
                                                 mat_T = res$estimates$mat_T,
                                                 means = res$estimates$means,
                                                 sdev = res$estimates$sd)
    # See how many states are recovered
    estimated_states <- sample_states_normal_(mat_R = post_prob,
                                              pi = vpi,
                                              mat_T = res$estimates$mat_T,
                                              means = res$estimates$means,
                                              sdev = res$estimates$sd,
                                              data = trial_data)
  }
  
  if (inherits(res, "hmm_mcmc_poisson")) {
    # Calculate posterior probabilities of hidden states
    post_prob <-  posterior_probabilities_poisgamma(data = trial_data,
                                                    pi = vpi,
                                                    mat_T = res$estimates$mat_T,
                                                    betas = res$estimates$betas,
                                                    alpha = res$estimates$alpha)
    # See how many states are recovered
    estimated_states <- sample_states_pois_(mat_R = post_prob,
                                            pi = vpi,
                                            mat_T = res$estimates$mat_T,
                                            betas = res$estimates$betas,
                                            alpha = res$estimates$alpha,
                                            data = trial_data)
  }
  
  m_multi <- list("target" = trial_states,
                  "prediction" = estimated_states)
  
  conf_mat <- cvms::confusion_matrix(targets = m_multi$target,
                                     predictions = m_multi$prediction)
  
  if (plot) {
    suppressWarnings(
      cvms::plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]],
                                  add_sums = TRUE)
    )
  }
  conf_mat
}


# hmm_simulate_normal_data_vary_K
hmm_simulate_normal_data_vary_K <- function(L, tr, sp, sigma0, maxK) {
  
  lmatT <- list()
  lvecMeans <- list()
  lvecEmit <- list()
  lvecStates <- list()
  
  for (k in 1:maxK) {
    vU <- rep(tr,k)
    vL <- rep(tr,k)
    lmatT[[k]] <- get_mat_T_(vU,vL)
    lvecMeans[[k]] <- seq(0, k, 1) / k * sp
    simData <- hmm_simulate_normal_data(L, lmatT[[k]], lvecMeans[[k]], sigma0)
    lvecEmit[[k]] <- simData$data
    lvecStates[[k]] <- simData$states
  }
  list("transition_matrices" = lmatT,
       "true_means" = lvecMeans,
       "simulated_data" = lvecEmit,
       "simulated_states" = lvecStates)
}


# hmm_mcmc_normal_vary_K
hmm_mcmc_normal_vary_K <- function(data,
                                   matT.init.list = FALSE,
                                   matT.sim.list,
                                   shape = 25,
                                   means.init.list = FALSE,
                                   sd.init = 2,
                                   iter = 500) {
  
  L <- length(data)
  
  if (is.vector(data) != TRUE) {
    stop("hmm_mcmc_normal_vary_K(): data input must be a vector", call. = FALSE)
  }
  
  maxK <- length(matT.sim.list)
  
  if (matT.init.list != FALSE) {
    maxK <- length(matT.init.list)
    if (length(means.init.list) != maxK) {
      stop("hmm_mcmc_normal_vary_K(): `matT.sim.list` and `means.init.list` must be the same length", call. = FALSE)
    }
  }
  
  if(length(sd.init) != 1 ) {
    stop("hmm_mcmc_normal_vary_K(): sd.init must be a constant", call. = FALSE)
  }
  
  ResK <- list()
  
  if (matT.init.list == FALSE) {
    for (l in 1:(maxK)) {
      matT.init <- matT.sim.list[[l]]
      
      for (i in 1:(l + 1)) {
        # matT.init[i,] <- MCMCprecision::rdirichlet(1, matT.init[i, ] * shape)
        matT.init[i,] <- gtools::rdirichlet(1, matT.init[i, ] * shape)
      }
      
      mn <- mean(data)
      std <- stats::sd(data) / (l + 1)
      means.init <- sort(stats::rnorm(l + 1, mean = mn, sd = std))
      hmm_mcmc_normal.out <- hmm_mcmc_normal(data = data,
                                             init_T = NULL,
                                             init_means = NULL,
                                             init_sd = NULL,
                                             prior_T = matT.init,
                                             prior_means = means.init,
                                             prior_sd = sd.init,
                                             iter = iter
      )
      ResK[[l]] <- hmm_mcmc_normal.out
    }
  } else {
    for (l in 1:(maxK - 1)) {
      hmm_mcmc_normal.out <- hmm_mcmc_normal(data = data,
                                             init_T = NULL,
                                             init_means = NULL,
                                             init_sd = NULL,
                                             prior_T = matT.init.list[[l]],
                                             prior_means = means.init.list[[l]],
                                             prior_sd = sd.init,
                                             iter = iter)
      
      ResK[[l]] <- hmm_mcmc_normal.out
    }
  }
  ResK
}


