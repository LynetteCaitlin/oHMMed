#' @references
#' Claus Vogl, Mariia Karapetiants, Burçin Yıldırım, Hrönn Kjartansdóttir, Carolin Kosiol, Juraj Bergman, Michal Majka, Lynette Caitlin Mikula.
#' Inference of genomic landscapes using ordered Hidden Markov Models with emission densities (oHMMed).
#' BMC Bioinformatics 25, 151 (2024). \doi{10.1186/s12859-024-05751-4}
#' @aliases oHMMed-package NULL

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
#' # Data stored in the object
#' plot(density(example_hmm_mcmc_normal$data), main = "")
#' 
#' # Priors used in simulation
#' example_hmm_mcmc_normal$priors
#' 
#' # Model
#' example_hmm_mcmc_normal
#' 
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
#' # Data stored in the object
#' hist(example_hmm_mcmc_gamma_poisson$data, 
#'      breaks = 50, xlab = "", main = "")
#' 
#' # Priors used in simulation
#' example_hmm_mcmc_gamma_poisson$priors
#' 
#' # Model
#' example_hmm_mcmc_gamma_poisson
#' 
#' summary(example_hmm_mcmc_gamma_poisson)
NULL