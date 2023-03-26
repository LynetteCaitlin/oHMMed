# Hidden Markov Model simulation with Poisson-gamma data


#' Calculate a Kullback-Leibler Divergence for a Discrete Distribution
#' 
#' @param p (numeric) probabilities
#'
#' @param q (numeric) probabilities
#'
#' @details
#' The Kullback-Leibler divergence for a discrete distribution
#' is calculated as follows:
#' 
#' \deqn{\sum_{i=1}^n p_i \log\Big(\frac{p_i}{q_i}\Big)}
#'
#' @return
#' Numeric vector
#'
#' @export
#'
#' @examples
#' # Simulate n Poisson distributed variates
#' n <- 1000
#' dist1 <- rpois(n, lambda = 1)
#' dist2 <- rpois(n, lambda = 5)
#' dist3 <- rpois(n, lambda = 20)
#' 
#' # Generate common factor levels
#' x_max <- max(c(dist1, dist2, dist3))
#' all_levels <- 0:x_max
#' 
#' # Estimate probability mass functions 
#' pmf_dist1 <- table(factor(dist1, levels = all_levels)) / n
#' pmf_dist2 <- table(factor(dist2, levels = all_levels)) / n
#' pmf_dist3 <- table(factor(dist3, levels = all_levels)) / n
#' 
#' # Visualise PMFs
#' barplot(pmf_dist1, col = "green", xlim = c(0, x_max))
#' barplot(pmf_dist2, col = "red", add = TRUE)
#' barplot(pmf_dist3, col = "blue", add = TRUE)
#' 
#' # Calculate distances
#' kullback_leibler_disc(pmf_dist1, pmf_dist2)
#' kullback_leibler_disc(pmf_dist1, pmf_dist3)
#' kullback_leibler_disc(pmf_dist2, pmf_dist3)

kullback_leibler_disc <- function(p, q) {
  
  if (!is.numeric(p)) {
    stop("kullback_leibler_disc(): `p` must be numeric", call. = FALSE)
  }
  if (!is.numeric(q)) {
    stop("kullback_leibler_disc(): `q` must be numeric", call. = FALSE)
  }
  if (any(p < 0) | any(p > 1)) {
    warning("kullback_leibler_disc(): `p` must be between 0 and 1", call. = FALSE)
  }
  if (any(q < 0) | any(q > 1)) {
    warning("kullback_leibler_disc(): `q` must be between 0 and 1", call. = FALSE)
  }
  
  if (length(p) != length(q)) {
    stop("kullback_leibler_disc(): `p` must be the same length as `q`", call. = FALSE)
  }
  
  p <- p + 1e-12
  q <- q + 1e-12
  
  sum(p * log(p / q))
}


#' Simulate Data Based on a Gamma-Poisson Model for a Hidden Markov Model Simulation
#'
#' @param L (integer) number of simulations
#'
#' @param mat_T (matrix) a square matrix with the initial state
#'
#' @param betas (numeric) \code{rate} parameter in \code{\link{rgamma}} for emission probabilities
#'
#' @param alpha (numeric) \code{shape} parameter in \code{\link{rgamma}} for emission probabilities
#'
#' @return
#' Returns a data vector "data", the "true" hidden states "states" used to generate the data vector
#' and prior probability of states "pi".
#'
#' @export
#'
#' @examples
#' mat_T <- rbind(c(1-0.01, 0.01, 0),
#'                c(0.01, 1-0.02, 0.01),
#'                c(0, 0.01, 1-0.01))
#' L <- 2^7
#' betas <- c(0.1, 0.3, 0.5)
#' alpha <- 1
#' 
#' sim_data <- hmm_simulate_gamma_poisson_data(L = L,
#'                                             mat_T = mat_T,
#'                                             betas = betas,
#'                                             alpha = alpha)
#' hist(sim_data$data, breaks = 40,
#'      main = "Histogram of Simulated Gamma-Poisson Data", xlab = "")
#' sim_data

hmm_simulate_gamma_poisson_data = function(L, mat_T, betas, alpha) {
  
  if (length(alpha) != 1) {
    stop("hmm_simulate_gamma_poisson_data(): shape parameter `alpha` must be of length 1", call. = FALSE)
  }
  
  if (alpha < 0) {
    stop("hmm_simulate_gamma_poisson_data(): shape parameter `alpha` must be non-negative", call. = FALSE)
  }
  
  if (any(betas < 0)) {
    stop("hmm_simulate_gamma_poisson_data(): rate parameters `betas` must be non-negative", call. = FALSE)
  }
  
  if (!is.matrix(mat_T)) {
    stop("hmm_simulate_gamma_poisson_data(): `mat_T` must be a numeric matrix", call. = FALSE)
  }
  
  if (any(is_row_sum_one_(mat_T) == FALSE)) {
    stop("hmm_simulate_gamma_poisson_data(): rows in the transition matrix `mat_T` must sum up to 1", call. = FALSE)
  }
  
  if (length(L) != 1) {
    stop("hmm_simulate_gamma_poisson_data(): the number of simulations `L` must be a single integer", call. = FALSE)
  }
  
  vstate <- rep(NA, L)
  vemit <- rep(NA, L)
  
  n_states <- nrow(mat_T)
  pi <- get_pi_(mat_T)
  
  vstate[1] <- sample.int(n = n_states, size = 1, replace = TRUE, prob = pi)
  vemit[1] <- stats::rpois(1, stats::rgamma(1, shape = alpha, rate = betas[vstate[1]]))
  
  for(i in 2:L) {
    p_state <- mat_T[vstate[i - 1], ]
    vstate[i] <- sample.int(n = n_states, size = 1, replace = TRUE, prob = p_state)
    vemit[i] <- stats::rpois(1, stats::rgamma(1, shape = alpha, rate = betas[vstate[i]]))
  }
  
  list("data" = vemit,
       "states" = factor(vstate, levels = 1:n_states),
       "pi" = pi)
}


# dNegBinom=function(y,b,a){
#     retval=-lfactorial(y)+lgamma(y+a)-lgamma(a)+a*log(b/(1+b))-y*log(1+b)
#     return(exp(retval))
# }
# <==>
# dnbinom(x = y, prob = b / (1 + b), size = a)


#' Forward-Backward Algorithm to Calculate the Posterior Probabilities of Hhidden States in Poisson-Gamma Model
#'
#' Forward-Backward Algorithm to Calculate the Posterior Probabilities of Hhidden States in Poisson-Gamma Model
#'
#' @param data (numeric) Poisson data
#'
#' @param pi (numeric) prior probability of states
#'
#' @param mat_T (matrix) transition probability matrix
#'
#' @param betas (numeric) vector with prior rates
#'
#' @param alpha (numeric) prior scale
#'
#' @details
#' TODO: Here details on how the calculation is made
#'
#' @references
#' TODO: Here some references
#'
#' @return
#' TODO: List with posterior probabilities (improve description)
#'
#' @export
#'
#' @examples
#' mat_T <- rbind(c(1-0.01,0.01,0),
#'                c(0.01,1-0.02,0.01),
#'                c(0,0.01,1-0.01))
#' L <- 2^10
#' betas <- c(0.1, 0.3, 0.5)
#' alpha <- 1
#'
#' sim_data <- hmm_simulate_gamma_poisson_data(L = L,
#'                                             mat_T = mat_T,
#'                                             betas = betas,
#'                                             alpha = alpha)
#' pi <- sim_data$pi
#' hmm_poison_data <- sim_data$data
#' hist(hmm_poison_data, breaks = 100)
#'
#' # Calculate posterior probabilities of hidden states
#' post_prob <- posterior_prob_gamma_poisson(data = hmm_poison_data,
#'                                           pi = pi,
#'                                           mat_T = mat_T,
#'                                           betas = betas,
#'                                           alpha = alpha)
#' str(post_prob)

posterior_prob_gamma_poisson <- function(data, pi, mat_T, betas, alpha) {
  
  cap_ <- .Machine$double.xmax
  floor_ <- .Machine$double.xmin
  
  L <- length(data)
  n_states <- nrow(mat_T)
  mat_F <- matrix(ncol = n_states, nrow = L)
  vs <- rep(1, L)
  
  pi[pi < floor_] <- floor_
  nb_prob <- betas / (1 + betas)
  
  vtemp <- cap_floor_(pi * stats::dnbinom(x = data[1], prob = nb_prob, size = alpha), cap_, floor_)
  vs[1] <- sum(vtemp)
  mat_F[1,] <- vtemp / vs[1]
  
  # Forward
  for (l in 2:L) {
    probs <- cap_floor_(stats::dnbinom(x = data[l], prob = nb_prob, size = alpha), cap_, floor_)
    vtemp <- mat_F[l - 1,] %*% mat_T * probs
    vs[l] <- sum(vtemp)
    mat_F[l, ] <- vtemp / vs[l]
  }
  
  # Backward
  mat_B <- matrix(nrow = L, ncol = n_states)
  mat_B[L, ] <- rep(1, n_states) / vs[L]
  
  for (l in L:2) {
    probs <- cap_floor_(stats::dnbinom(x = data[l], prob = nb_prob, size = alpha), cap_, floor_)
    mat_B[l - 1, ] <- mat_T %*% (probs * mat_B[l, ]) / vs[l - 1]
  }
  
  list(F = mat_F, B = mat_B, s = vs)
}


#' @keywords internal
sample_states_pois_ <- function(mat_R, mat_T, pi, betas, alpha, data) {
  
  L <- length(mat_R$F[ ,1])
  n_states <- length(mat_R$F[1, ])
  states <- rep(NA, L)
  cap_ <- .Machine$double.xmax
  floor_ <- .Machine$double.xmin
  
  p <- cap_floor_(mat_R$F[1, ] * mat_R$F[1, ] * mat_R$s[1], cap_, floor_)
  p <- p / sum(p)
  p[is.na(p)] <- floor_
  # states[1] <- sample(1:n_states, size = 1, prob = p) # will exactly reproduce old function, 3 lines below must be commented out.
  all_runif <- stats::runif(L) # generating L unifs at once is faster than generating L times 1 unif
  states[1] <- which.max(cumsum(p) - all_runif[1] > 0) # new faster version
  # states[1] <- which.max(cumsum(p) - stats::runif(1) > 0) # old version
  
  
  nb_prob <- betas / (1 + betas)
  probs <- matrix(NA, nrow = L, ncol = length(betas))
  for (i in 1:L) {
    probs[i, ] <- nb_prob^alpha * (1 - nb_prob)^data[i]
  }
  probs <- cap_floor_(probs, cap_, floor_)
  
  for (l in 2:L) {
    p <- mat_T[states[l - 1], ] * probs[l, ] * mat_R$B[l, ]
    p <- p / sum(p)
    p[is.na(p)] <- floor_
    # states[l] <- sample(1:n_states, size = 1, prob = p) # will exactly reproduce old function
    states[l] <- which.max(cumsum(p) - all_runif[l] > 0) # new faster version
    # states[l] <- which.max(cumsum(p) - stats::runif(1) > 0) # old version
  }
  factor(states, levels = 1:n_states)
}


#' @keywords internal
sample_betas_alpha_ <- function(cur_betas, cur_alpha, prior_betas, prior_alpha, states = NULL, data = NULL) {
  
  N <- length(data)
  betas <- NULL
  alpha <- NULL
  n_states <- length(cur_betas)
  loglambda_p <- 0 # note: loglambda_p = log(lambda_p) in script
  
  #### begin NEW LCM & CV: no data -> sample betas from priors and alphas from close to prior
  if (is.null(states) & is.null(data)) {
    
    # betas <- sort(cur_betas, decreasing = TRUE) # MM: BUG (when applying "init_hmm_mcmc_gamma_poisson_" the cur_betas=NULL --> sorts(NULL) --> NULL. No need for sorting)
    # betas[i] <- stats::rgamma(1, shape = 1, rate = 1 / prior_betas) # MM: BUG (in addition to the above: missing for-loop. 3 betas with different rate parameters can be simulated in a vectorized way though)
    betas <- stats::rgamma(length(prior_betas), shape = 1, rate = 1 / prior_betas)
    
    # shape_param <- prior_alpha
    alpha <- stats::rgamma(1, shape = prior_alpha, rate = 1)
    #### end NEW LCM & CV
    
  } else {
    
    lvl <- list()
    for (i in 1:n_states) {
      vl <- data[states == i]
      N_i <- length(vl)
      vl <- stats::rgamma(N_i, shape = cur_alpha + vl, rate = cur_betas[i] + 1)
      lvl[[i]] <- vl
      rate_gamma <- 1
      post_shape <- cur_alpha * N_i + 1
      post_rate <- sum(lvl[[i]]) + 1 / prior_betas[i]
      ##print(c(post_shape, post_rate, post_shape / post_rate))
      betas[i] <- stats::rgamma(1, shape = post_shape, rate = post_rate)
    }
    
    betas <- sort(betas, decreasing = TRUE) 
    for (i in 1:n_states) {
      if (N > 0) {
        loglambda_p <- loglambda_p + sum(log(betas[i] * lvl[[i]]))
      }
    }  
    
    HM <- exp(loglambda_p / (N + prior_alpha))
    shape_param <- (N + prior_alpha) * cur_alpha
    
    alpha_star <- stats::rgamma(1, shape = shape_param, rate = N)
    lh_star <- 0
    lh_old <- 0
    for (i in 1:n_states) {
      lh_star <- lh_star + (log(HM) * alpha_star - lgamma(alpha_star)) * N
      lh_old <- lh_old + (log(HM) * cur_alpha - lgamma(cur_alpha)) * N
    }
    
    d1 <- stats::dgamma(alpha_star, shape = shape_param, rate = N, log = TRUE)
    d2 <- stats::dgamma(cur_alpha, shape = shape_param, rate = N, log = TRUE)
    lr <- lh_star - d1 - (lh_old - d2)
    alpha <- cur_alpha
    
    ## print(c(HM, N, alpha_star, cur_alpha, d1, d2, lh_star, lh_old, lr))
    
    if (stats::runif(1) <= exp(lr)) {
      alpha <- alpha_star
    }
  }
  ##print(c(betas, alpha))
  list(betas = betas, alpha = alpha)
}


#' @keywords internal
init_hmm_mcmc_gamma_poisson_ <- function(data, prior_T, prior_betas, prior_alpha,
                                         init_T, init_betas, init_alpha, 
                                         verbose, iter, warmup, thin, chain_id = NULL) {
  
  if (iter < 2) {
    stop("hmm_mcmc_gamma_poisson(): `iter` needs to be bigger than 1", call. = FALSE)
  }
  
  if (!is.integer(data)) {
    stop("hmm_mcmc_gamma_poisson(): `data` needs to be an integer vector", call. = FALSE)
  }
  
  if (sum(is.na(data)) > 0) {
    stop("hmm_mcmc_gamma_poisson(): `data` contains missing values", call. = FALSE)
  }
  
  if (any(is_row_sum_one_(prior_T) == FALSE)) {
    stop("hmm_mcmc_gamma_poisson(): rows in the transition matrix `prior_T` must sum up to 1", call. = FALSE)
  }
  if (length(prior_betas) != nrow(prior_T) | length(prior_betas) != ncol(prior_T)) {
    stop("hmm_mcmc_gamma_poisson(): number of states is not the same between input variables", call. = FALSE)
  }
  
  # the algorithm takes care of it anyway....
  #is_prior_beta_decreasing <- all(sort(prior_betas, decreasing = TRUE) == prior_betas)
  #if (!is_prior_beta_decreasing) {
  #  warning("hmm_mcmc_gamma_poisson(): `prior_betas` should be sorted in decreasing order", call. = FALSE)
  #}
  
  if (length(prior_alpha) != 1) {
    stop("hmm_mcmc_gamma_poisson(): `prior_alpha` must be of length 1", call. = FALSE)
  }
  if (warmup >= iter + 2) {
    stop("hmm_mcmc_gamma_poisson(): `warmup` must be lower than `iter`", call. = FALSE)
  }
  if (thin > iter - warmup) {
    stop("hmm_mcmc_gamma_poisson(): `thin` cannot exceed iterations after warmup period", call. = FALSE)
  }
  
  # Initialization
  out <- sample_betas_alpha_(cur_betas = init_betas, cur_alpha = init_alpha,  
                             prior_betas = prior_betas, prior_alpha = prior_alpha) ## CV: currently, pior_alpha and prior_betas useless
  if (is.null(init_betas)) {
    init_betas <- out$betas
  }
  
  if (is.null(init_alpha)) {
    init_alpha <- out$alpha
  }
  
  if (is.null(init_T)) {
    init_T <- sample_T_(prior_mat = prior_T)
  }
  
  if (length(init_betas) != nrow(init_T) | length(init_betas) != ncol(init_T)) {
    stop("hmm_mcmc_gamma_poisson(): number of states is not the same between input variables", call. = FALSE)
  }
  
  if (length(init_alpha) != 1) {
    stop("hmm_mcmc_gamma_poisson(): `init_alpha` must be of length 1", call. = FALSE)
  }
  
  if (any(is_row_sum_one_(init_T) == FALSE)) {
    stop("hmm_mcmc_gamma_poisson(): rows in the transition matrix `init_T` must sum up to 1", call. = FALSE)
  }
  # the algorithm takes care of it anyway....
  #is_init_betas_decreasing <- all(sort(init_betas, decreasing = TRUE) == init_betas)
  #if (!is_init_betas_decreasing) {
  #  warning("hmm_mcmc_gamma_poisson(): `init_betas` should be sorted in decreasing order", call. = FALSE)
  #}
  
  #NEW LCM
  lambda1 <- sum(data) / length(data)
  lambda2 <- numeric(length(init_betas)) # rep(NA, length(init_betas))
  for (i in 1:length(init_betas)) { 
    lambda2[i] <- mean(stats::rgamma(1000, init_alpha, init_betas[i]))
  }
  range1 <- lambda1 < lambda2
  range2 <- lambda1 > lambda2
  # END NEW LCM
  # TODO: estimate more or less optimal range for prior betas
  # sample_mean <- mean(data)
  # sample_sd <- stats::sd(data)
  # abs_mean_ratios <- abs(init_means) / abs(sample_mean)
  # sd_ratio <- init_sd / sample_sd

  chain_char <- NULL
  if (!is.null(chain_id)) {
    chain_char <- paste0("(chain ", chain_id, ") ")
  }
  
  if (verbose) {
     set_alpha <- (mean(data)^2) / ((stats::var(data) - mean(data)) / length(init_betas))
     if( (init_alpha > (set_alpha + 3)) | (init_alpha < (set_alpha - 3)) ){
       message("hmm_mcmc_gamma_poisson(): ", chain_char, "initial alpha is not close to the recommended value")
     }
     
    if ((sum(range1) == length(lambda2)) | (sum(range2) == length(lambda2))) {
      message("hmm_mcmc_gamma_poisson(): ", chain_char, "rate parameter of observed distribution is either above or below all of the initial rate parameters")
    }
    if (any(lambda2 > stats::var(data))) {
      message("hmm_mcmc_gamma_poisson(): ", chain_char, "at least one initial rate parameter is greater than the overall variance")
    }
  }
  
  init_pi <- get_pi_(init_T)
  
  init_mat_res <- posterior_prob_gamma_poisson(data = data,
                                               pi = init_pi,
                                               mat_T = init_T,
                                               betas = init_betas,
                                               alpha = init_alpha)
  
  init_states <- sample_states_pois_(mat_R = init_mat_res,
                                     pi = init_pi,
                                     mat_T = init_T,
                                     betas = init_betas,
                                     alpha = init_alpha,
                                     data = data)
  
  list(init_states = init_states, init_betas = init_betas,
       init_alpha = init_alpha, init_T = init_T, init_pi = init_pi,
       init_mat_res = init_mat_res)
}


#' MCMC Simulation of a Hidden Markov Gamma-Poisson Model
#'
#' @param data (numeric) data
#'
#' @param prior_T (matrix) prior transition matrix
#'
#' @param prior_betas (numeric) prior beta parameters
#'
#' @param prior_alpha (numeric) a single prior alpha parameter. By default, \code{prior_alpha=1}
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
#' @param init_betas (numeric) \code{optional parameter}; initial beta parameters
#'
#' @param init_alpha (numeric) \code{optional parameter}; initial alpha parameter
#'
#' @param print_params (logical) \code{optional parameter}; print estimated parameters every iteration. By default, \code{TRUE}
#'
#' @param verbose (logical) \code{optional parameter}; print additional messages. By default, \code{TRUE}
#'
#' @details
#' TODO: Here details
#'
#' @references
#' TODO: Here references
#'
#' @return
#' List with following elements:
#' \itemize{
#'   \item data: data used for simulation
#'   \item samples: list with samples
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
#' # Simulate Poisson-Gamma data
#' N <- 2^10
#' true_T <- rbind(c(0.95, 0.05, 0),
#'                 c(0.025, 0.95, 0.025),
#'                 c(0.0, 0.05, 0.95))
#'
#' true_betas <- c(2, 1, 0.1)
#' true_alpha <- 1
#'
#' simdata_full <- hmm_simulate_gamma_poisson_data(L = N,
#'                                                 mat_T = true_T,
#'                                                 betas = true_betas,
#'                                                 alpha = true_alpha)
#' simdata <- simdata_full$data
#' hist(simdata, breaks = 40, probability = TRUE,  
#'      main = "Distribution of the simulated Poisson-Gamma data")
#' lines(density(simdata), col = "red")
#' 
#' # Set numbers of states to be inferred
#' n_states_inferred <- 3
#' 
#' # Set priors
#' prior_T <- generate_random_T(n_states_inferred)
#' prior_betas <- c(1, 0.5, 0.1)
#' prior_alpha <- 3
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
#' res <- hmm_mcmc_gamma_poisson(data = simdata,
#'                               prior_T = prior_T,
#'                               prior_betas = prior_betas,
#'                               prior_alpha = prior_alpha,
#'                               iter = iter,
#'                               warmup = warmup,  
#'                               thin = thin,
#'                               seed = seed,
#'                               print_params = print_params,
#'                               verbose = verbose)
#' res
#'
#' summary(res)# summary output can be also assigned to a variable
#' 
#' coef(res) # extract model estimates
#' 
#' # plot(res) # MCMC diagnostics

hmm_mcmc_gamma_poisson <- function(data,
                                   prior_T,
                                   prior_betas,
                                   prior_alpha = 1,
                                   iter = 1500,
                                   warmup = floor(iter / 1.5),
                                   thin = 1,
                                   seed = sample.int(.Machine$integer.max, 1),
                                   init_T = NULL,
                                   init_betas = NULL,
                                   init_alpha = NULL,
                                   print_params = TRUE,
                                   verbose = TRUE) {
  
  set.seed(seed)
  
  init_data <- init_hmm_mcmc_gamma_poisson_(data, prior_T, prior_betas, prior_alpha,
                                            init_T, init_betas, init_alpha,
                                            verbose, iter, warmup, thin)
  
  n_data <- length(data)
  n_states <- length(prior_betas)
  
  mean_states <- matrix(0, ncol = n_states, nrow = n_data)
  all_betas <- matrix(0, ncol = n_states, nrow = iter)
  all_alpha <- numeric(iter)
  all_mat_T <- array(dim = c(dim(prior_T), iter))
  vlh <- numeric(iter)
  
  all_betas[1, ] <- init_data$init_betas
  all_alpha[1] <- init_data$init_alpha
  all_mat_T[, ,1] <- init_data$init_T
  vlh[1] <- sum(log(init_data$init_mat_res$s))
  states <- init_data$init_states
  betas <- init_data$init_betas
  alpha <- init_data$init_alpha
  
  # Run sampler
  for (it in 2:iter) {
    
    print_progress_(it, iter, verbose = verbose)
    m <- sample_betas_alpha_(cur_betas = betas,
                             cur_alpha = alpha,
                             prior_betas = prior_betas,
                             prior_alpha = prior_alpha,
                             states = states,
                             data = data)
    
    betas <- m$betas
    alpha <- m$alpha
    mat_T <- sample_T_(states, prior_mat = prior_T)
    pi <- get_pi_(mat_T)
    mat_res <- posterior_prob_gamma_poisson(data = data,
                                            pi = pi,
                                            mat_T = mat_T,
                                            betas = betas,
                                            alpha = alpha)
    
    states <- sample_states_pois_(mat_R = mat_res,
                                  pi = pi,
                                  mat_T = mat_T,
                                  betas = betas,
                                  alpha = alpha,
                                  data = data)
    if (iter > warmup) {
      for (i in 1:n_states) {
        mean_states[ ,i] <- mean_states[ ,i] + (states == i)
      }
    }
    
    all_betas[it, ] <- betas
    all_alpha[it] <- alpha
    all_mat_T[ , ,it] <- mat_T
    vlh[it] <- sum(log(mat_res$s))
    
    if (print_params) {
      print(c(colMeans(all_betas[1:it, ]), mean(all_alpha[1:it])))
    }
  }
  
  # Prepare outputs
  idx <- seq.int(warmup + 1, by = thin, to = iter)
  
  colnames(all_betas) <- paste0("beta[", 1:n_states,"]")
  all_means <- all_alpha / all_betas
  colnames(all_means) <- paste0("means[", 1:n_states,"]")
  
  samples <- list(betas = all_betas,
                  alpha = all_alpha,
                  means = all_means,
                  mat_T = all_mat_T,
                  vlh = vlh)
  
  posterior_states <- factor(max.col(mean_states, "first"), 1:n_states)
  
  estimates <- list(betas = colMeans(all_betas[idx, ]),
                    alpha = mean(all_alpha[idx]),
                    means = colMeans(all_means[idx, ]),
                    mat_T = apply(all_mat_T[ , ,idx], c(1,2), mean),
                    posterior_states = posterior_states,
                    posterior_states_prob = mean_states / iter,
                    log_likelihood = vlh[idx])    
  
  priors <- list(prior_betas = prior_betas,
                 prior_alpha = prior_alpha,
                 prior_T = prior_T)
  
  inits <- list(init_states = init_data$init_states,
                init_betas = init_data$init_betas,
                init_alpha = init_data$init_alpha,
                init_T = init_data$init_T)
  
  last_iter <- list(betas = all_betas[iter, ],
                    alpha = all_alpha[iter],
                    means = all_means[iter, ],
                    mat_T = all_mat_T[ , ,iter])
  
  info <- list(model_name = "hmm_mcmc_gamma_poisson",
               date = as.character(Sys.time()),
               seed = seed,
               iter = iter,
               warmup = warmup,
               thin = thin,
               n_states = length(prior_betas))
  
  res <- list(data = data,
              samples = samples,
              estimates = estimates,
              idx = idx,
              priors = priors,
              inits = inits,
              last_iter = last_iter,
              info = info)
  
  class(res) <- "hmm_mcmc_gamma_poisson"
  res$info$object_size <- format(utils::object.size(res), "Mb")
  res
}


#' @keywords internal
#' @export

print.hmm_mcmc_gamma_poisson <- function(x, ...) {
  info <- x$info
  cat("Model:", "HMM Gamma-Poisson", "\n")
  cat("Type:", "MCMC", "\n")
  cat("Iter:", info$iter, "\n")
  cat("Warmup:", info$warmup, "\n")
  cat("Thin:", info$thin, "\n")
  cat("States:", length(x$priors$prior_betas), "\n")
}


#' @keywords internal
#' @export

summary.hmm_mcmc_gamma_poisson <- function(object, ...) {
  
  info <- object$info
  idx <- object$idx
  data <- object$data
  
  beta_est <- object$estimates$betas
  alpha_est <- object$estimates$alpha
  T_est <- object$estimates$mat_T
  
  means_est <- object$estimates$means
  
  post_states <- object$estimates$posterior_states
  state_tab <- table(post_states, dnn = "")
  
  kl_list <- rep(NA, 500)
  for (j in 1:500) {
    sim_output <- unlist(lapply(1:length(state_tab), function(i) {
      stats::rnbinom(state_tab[i],
                     size = alpha_est,
                     prob = beta_est[i] / (1 + beta_est[i]))
    }))
    jth_levels <- 0:max(c(data, sim_output))
    dens_data_table <- table(factor(data, levels = jth_levels))
    dens_sim_table <- table(factor(sim_output, levels = jth_levels))
    
    dens_data <- dens_data_table / sum(dens_data_table)
    dens_sim <- dens_sim_table / sum(dens_sim_table)
    
    kl_list[j] <- kullback_leibler_disc(dens_data, dens_sim)
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
  #names(post_info) <- c("mean", "sd", "median")
  
  group_comparison <- rep(NA, length(beta_est) - 1)
  for (k in 1:(length(beta_est) - 1)) {
    gr1 <- object$data[object$estimates$posterior_states == k]
    gr2 <- object$data[object$estimates$posterior_states == (k + 1)]
    group_comparison[k] <- stats::poisson.test(x = c(sum(gr1), sum(gr2)), 
                                               T = c(length(gr1), length(gr2)), 
                                               alternative = "l")$p.value
  }
  
  summary_res <- list("estimated_betas" = beta_est,
                      "estimated_alpha" = alpha_est,
                      "estimated_means" = means_est,
                      "estimated_transition_rates" = T_est,
                      "assigned_states" = state_tab,
                      "approximate_kullback_leibler_divergence" = kl_div,
                      "log_likelihood" = ll_info,
                      "state_differences_significance" = group_comparison)
  
  #### NEW LCM: remove below????
  #"log_posterior" = post_info)
  
  cat("Estimated betas:\n")
  print(summary_res$estimated_betas)
  cat("\n")
  
  if (stats::sd(object$samples$alpha) > 0) {
    cat("Estimated alpha:\n")
    cat(summary_res$estimated_alpha)
    cat("\n")
    cat("\n")
  }
  
  cat("Estimated means:\n")
  cat(summary_res$estimated_means)
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
  
  cat("Significance of Difference between Rates (stepwise):\n")
  print(summary_res$state_differences_significance)
  cat("\n")
  
  #### NEW LCM: remove below????
  #cat("Log Posterior:\n")
  #print(summary_res$log_posterior)
  
  invisible(summary_res)
}


#' @rdname coef.hmm_mcmc_normal
#'
#' @export
#' @export coef.hmm_mcmc_gamma_poisson
#'
#' @examples
#' coef(example_hmm_mcmc_gamma_poisson)

coef.hmm_mcmc_gamma_poisson <- function(object, ...) {
  est <- object$estimates
  list(betas = est$betas,
       alpha = est$alpha,
       means = est$means,
       mat_T = est$mat_T)
}


#' Plot Diagnostics for \code{hmm_mcmc_gamma_poisson} Objects
#' 
#' This function creates a variety of diagnostic plots that can be useful when 
#' conducting Markov Chain Monte Carlo (MCMC) simulation of a gamma-poisson hidden Markov model (HMM). 
#' These plots will help to assess convergence, fit, and performance of the MCMC simulation
#'
#' @param x (hmm_mcmc_gamma_poisson) HMM MCMC gamma-poisson object
#'
#' @param simulation (logical) TODO:
#'
#' @param true_betas (numeric) true betas. To be used if \code{simulation=TRUE}
#'
#' @param true_alpha (numeric) true alpha. To be used if \code{simulation=TRUE}
#'
#' @param true_mat_T (matrix) \code{optional parameter}; true transition matrix. To be used if \code{simulation=TRUE}
#'
#' @param true_states (integer) \code{optional parameter}; true states. To be used if \code{simulation=TRUE}
#' 
#' @param show_titles (logical) if \code{TRUE} then titles are shown for all graphs. By default, \code{TRUE}
#'
#' @param ... not used
#'
#' @return
#' Several diagnostic plots that can be used to evaluate the MCMC simulation
#' of the gamma-poisson HMM
#'
#' @export
#' @export plot.hmm_mcmc_gamma_poisson
#'
#' @examples
#' \donttest{
#' plot(example_hmm_mcmc_gamma_poisson)
#' }

plot.hmm_mcmc_gamma_poisson <- function(x,
                                        simulation = FALSE,
                                        true_betas = NULL,
                                        true_alpha = NULL,
                                        true_mat_T = NULL,
                                        true_states = NULL,
                                        show_titles = TRUE,
                                        ...) {
  
  info <- x$info
  data <- x$data
  idx <- x$idx
  
  if (simulation) {
    cond <- any(c(is.null(true_betas), is.null(true_alpha), is.null(true_mat_T),
                  is.null(true_states)))
    if (cond) {
      stop("plot.hmm_mcmc_gamma_poisson(): if `simulation=TRUE` then `true_betas` ",
           "`true_alpha`, `true_mat_T` and `true_states` must be defined", call. = FALSE)
    }
  }
  
  # Diagnostics betas
  all_betas <- convert_to_ggmcmc(x, pattern = "beta")
  n_betas <- attributes(all_betas)$nParameters
  facet_means <- ggplot2::facet_wrap(~ Parameter, ncol = floor(n_betas / 2), scales = "free")
  
  mdens <- ggmcmc::ggs_density(all_betas) + 
    facet_means + 
    ggplot2::labs(x = "Value", 
                  y = "Density", 
                  title = if (show_titles) "Densities of Betas" else NULL)
  
  mtrace <- ggmcmc::ggs_traceplot(all_betas) + 
    facet_means + 
    ggplot2::labs(x = "Iteration", 
                  y = "Value", 
                  title = if (show_titles) "Traceplots of Betas" else NULL)
  
  # Diagnostics transitions
  all_T <- convert_to_ggmcmc(x, pattern = "T")
  n_t <- attributes(all_T)$nParameter
  facet_t <- ggplot2::facet_wrap(~ Parameter, ncol = sqrt(n_t), scales = "free")
  labels_t <- ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = "."))
  
  Ttrace <- ggmcmc::ggs_traceplot(all_T) + 
    facet_t + 
    labels_t + 
    ggplot2::labs(x = "Iteration", 
                  y = "Value",
                  title = if (show_titles) "Traceplots of Parameters in Transition Matrix" else NULL)
  
  Tdens <- ggmcmc::ggs_density(all_T) + 
    facet_t + 
    labels_t + 
    ggplot2::labs(x = "Value", 
                  y = "Density",
                  title = if (show_titles) "Densities of Parameters in Transition Matrix" else NULL)
  
  # Diagnostics alpha
  df_alpha <- convert_to_ggmcmc(x, "alpha")
  sdtrace <- ggmcmc::ggs_traceplot(df_alpha) + 
    ggplot2::labs(x = "Iteration", 
                  y = "Value", 
                  title = if (show_titles) "Traceplot of Alpha" else NULL)
  
  sddens <- ggmcmc::ggs_density(df_alpha) +
    ggplot2::labs(x = "Value", 
                  y = "Density",
                  title = if (show_titles) "Density of Alpha" else NULL)
  
  # Diagnostics mean
  all_means <- convert_to_ggmcmc(x, pattern = "means")
  n_means <- attributes(all_means)$nParameters
  facet_me <- ggplot2::facet_wrap(~ Parameter, ncol = floor(n_means / 2), scales = "free")
  
  medens <- ggmcmc::ggs_density(all_means) + 
    facet_me + 
    ggplot2::labs(x = "Value", 
                  y = "Density",
                  title = if (show_titles) "Densities of Means" else NULL)
  
  metrace <- ggmcmc::ggs_traceplot(all_means) + 
    facet_me + 
    ggplot2::labs(x = "Iteration", 
                  y = "Value",
                  title = if (show_titles) "Traceplots of Means" else NULL)
  
  # Likelihood trace
  lltrace <- x$estimates$log_likelihood
  lltrace_df <- as.data.frame(cbind(c((info$warmup + 1):info$iter), lltrace))
  names(lltrace_df) <- c("iteration", "log_likelihood")
  
  llplot <- ggplot2::ggplot(lltrace_df, ggplot2::aes_string(x = "iteration", 
                                                            y = "log_likelihood")) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iteration", 
                  y = "Log-likelihood",
                  title = if (show_titles) "Traceplot of Log-Likelihood" else NULL)
  
  # Confusion matrix
  if (simulation) {
    m_multi <- list("target" = true_states,
                    "prediction" = x$estimates$posterior_states)
    conf_mat <- cvms::confusion_matrix(targets = m_multi$target,
                                       predictions = m_multi$prediction)
    
    conf_mat_plot <- suppressWarnings(
      cvms::plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]],
                                  add_sums = TRUE) +
        ggplot2:::labs(title = if (show_titles) "Confusion Matrix" else NULL)
    )
  }
  
  # # Check assignment of states along chromosome
  states_df <- as.data.frame(cbind(1:length(data), data, x$estimates$posterior_states))
  post_means <- numeric(length(data))
  
  for (l in 1:length(data)) {
    post_means[l] <- sum(x$estimates$alpha / x$estimates$betas * x$estimates$posterior_states_prob[l, ])
  }
  
  states_df$post_means <- post_means
  names(states_df) <- c("position", "data", "posterior_states", "posterior_means")
  states_df$posterior_states <- as.factor(states_df$posterior_states)
  
  statesplot <- ggplot2::ggplot(states_df, ggplot2::aes_string(x = "position", y = "data")) +
    ggplot2::geom_line(col = "grey") +
    ggplot2::geom_point(ggplot2::aes_string(colour = "posterior_states"), shape = 20, size = 1.5, alpha = 0.75) +
    ggplot2::geom_line(ggplot2::aes_string(x = "position", y = "posterior_means"), size = 0.15) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Post States")) +
    ggplot2::labs(x = "Position", 
                  y = "Data",
                  title = if (show_titles) "States Plot" else NULL)
  
  if (simulation) {
    states_df2 <- as.data.frame(cbind(1:length(data), data, true_states))
    post_means <- (true_alpha / true_betas)[true_states]
    states_df2$post_means <- post_means
    names(states_df2) <- c("position", "data", "true_states", "true_means")
    states_df2$true_states <- as.factor(states_df2$true_states)
    
    statesplot2 <- ggplot2::ggplot(states_df2, ggplot2::aes_string(x = "position", y = "data")) +
      ggplot2::geom_line(col = "grey") +
      ggplot2::geom_point(ggplot2::aes_string(colour = "true_states"), shape = 20, size = 1.5, alpha = 0.75) +
      ggplot2::geom_line(ggplot2::aes_string(x = "position", y = "true_means"), size = 0.15) +
      ggplot2::guides(colour = ggplot2::guide_legend(title = "True States")) +
      ggplot2::labs(x = "Position", y = "Data")
  }
  
  # Check fit
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
    sim_output1 <- c(sim_output1, sim_output)
  }
  sim_output1 <- sim_output1[2:length(sim_output1)]
  dens_data1 <- table(factor(data, levels = 0:max(c(data,sim_output1)))) / sum(table(factor(data, levels = 0:max(c(data, sim_output1))))) 
  dens_sim1 <- table(factor(sim_output1, levels = 0:max(c(data, sim_output1)))) / sum(table(factor(sim_output1, levels = 0:max(c(data, sim_output1)))))
  
  kl_plot_f <- function() {
    vcd::rootogram(as.numeric(dens_data1), 
                   as.numeric(dens_sim1),
                   lines_gp = grid::gpar(col = "red", lwd = 1), 
                   main = if (show_titles) "Model Fit" else NULL,
                   points_gp = grid::gpar(col = "black"), 
                   pch = "",
                   ylab = "Frequency (sqrt)",
                   xlab = "Number of Occurrences")
  }
  
  dens_df <- as.data.frame(cbind(rep("observed", length(data)), data))
  
  names(dens_df) <- c("data_type", "value")
  dens_df$value <- as.numeric(dens_df$value)
  
  klplot <- ggplot2::ggplot(dens_df, ggplot2::aes_string(x = "value")) +
    ggplot2::geom_histogram(bins = floor(dim(table(factor(data, levels = 0:max(data))))), 
                            fill = "grey", position = "identity") +
    ggplot2::scale_color_manual(values = c("black")) +
    ggplot2::geom_vline(xintercept = x$estimates$means, color = "black", size = 0.2) +
    ggplot2::labs(x = "Number of Occurences", 
                  y = "Frequency",
                  title = if (show_titles) "Observed Counts and Inferred Means" else NULL) 
  
  if (simulation) {
    klplot <- ggplot2::ggplot(dens_df, ggplot2::aes_string(x = "value")) +
      ggplot2::geom_histogram(bins = floor(dim(table(factor(data, levels = 0:max(data))))), 
                              fill = "grey", position = "identity") +
      ggplot2::scale_color_manual(values = c("black")) +
      ggplot2::geom_vline(xintercept = x$estimates$means, color = "black", size = 0.6) +
      ggplot2::geom_vline(xintercept = true_alpha /true_betas, color = "blue", 
                          size = 0.4, linetype = "dotted") + 
      ggplot2::labs(x = "Number of Occurences", 
                    y = "Frequency",
                    title = if (show_titles) "Observed Counts and Inferred (and True) Means" else NULL)
    
  }
  
  if (simulation) {
    plotlist <- list(mtrace, mdens, Ttrace, Tdens, sdtrace, sddens, metrace, 
                     medens, llplot, conf_mat_plot, klplot)
  } 
  
  if (isFALSE(simulation)) {
    plotlist <- list(mtrace, mdens, Ttrace, Tdens, sdtrace, sddens, metrace, 
                     medens, llplot, statesplot, klplot)
  } 
  
  for (ii in 1:length(plotlist)) {
    
    if (simulation && ii == (length(plotlist) - 2)) {
      gridExtra::grid.arrange(statesplot, statesplot2, nrow = 2)
    }
    invisible(utils::capture.output(print(plotlist[[ii]])))
  }
  kl_plot_f()   
}
