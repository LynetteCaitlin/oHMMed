#' @keywords internal 
#' 
#' @references
#' TODO: Reference to the paper 
#' 
"_PACKAGE"


#' Example of a Simulated Normal Model
#'
#' @docType data
#' 
#' @keywords datasets
#' 
#' @name example_hmm_mcmc_normal
#' 
#' @usage example_hmm_mcmc_normal
#' 
#' @format hmm_mcmc_normal object
#' 
#' @examples
#' # Data
#' plot(density(example_hmm_mcmc_normal$data), main = "")
#' 
#' # Priors
#' example_hmm_mcmc_normal$priors
#' 
#' # Model
#' example_hmm_mcmc_normal
#' summary(example_hmm_mcmc_normal)
NULL


#' Example of a Simulated Gamma-Poisson Model
#'
#' @docType data
#' 
#' @keywords datasets
#' 
#' @name example_hmm_mcmc_gamma_poisson
#' 
#' @usage example_hmm_mcmc_gamma_poisson
#' 
#' @format hmm_mcmc_gamma_poisson object
#' 
#' @examples
#' # Data
#' hist(example_hmm_mcmc_gamma_poisson$data, 
#'      breaks = 50, xlab = "", main = "")
#' 
#' # Priors
#' example_hmm_mcmc_gamma_poisson$priors
#' 
#' # Model
#' example_hmm_mcmc_gamma_poisson
#' summary(example_hmm_mcmc_gamma_poisson)
NULL