#' @references
#' Inference of Genomic Landscapes using Ordered Hidden Markov Models with Emission Densities (oHMMed)
#' Claus Vogl, Mariia Karapetiants, Burçin Yildirim, Hrönn Kjartansdottir, Carolin Kosiol, Juraj Bergman, Michal Majka, Lynette Caitlin Mikula
#' bioRxiv 2023.06.26.546495; doi: \url{https://doi.org/10.1101/2023.06.26.546495}
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