#' Compute likelihood weights
#'
#' This function computes the \emph{normalized} state weights and the \emph{unnormalized}
#' parameter weight for a given parameter vector and its corresponding states particles.
#'
#' @param states A data frame with columns \code{psi}, \code{I} and \code{R}
#' (corresponding to \eqn{\psi, i} and \eqn{r}, respectively).
#' @param fixed_inputs Necessary arguments to further specify the chosen
#' transition kernel, a list. Given the updating schemes describes in the paper,
#' it should include the following:
#' \itemize{
#' \item \code{tot_pop}: total number of individuals in the considered population
#' (indicated by \eqn{N} in the paper), a positive integer;
#' \item \code{g}: inverse of the average time a non-detected individual stays
#' infectious (indicated by \eqn{\gamma} in the paper), a positive scalar;
#' \item \code{q}: probability of detecting an infectious person through testing
#' (indicated by \eqn{q} in the paper), a non-negative scalar;
#' \item \code{d_f}: inverse of the average time an infected person enjoys immunity
#' infectious (indicated by \eqn{\delta} in the paper), a positive scalar;
#' \item \code{quarantine}: logical; if TRUE, the transition kernel corresponding
#' to the model with quarantine is used.
#' }
#' Defaults to \code{stochSIR_state_fixed_inputs}, which is available to users upon
#' package loading.
#' @param current_obs Number of people detected at the current time point, a
#' positive integer.
#' @return A list containing \code{state_norm_weights} (normalized state weights)
#' and \code{param_weight} (unnormalized parameter weight for the parameter vector
#' corresponding to the \code{states} data frame passed as argument).
#' @export
stochSIR_likelihood_weights <- function(states, current_obs, fixed_inputs){
  log_state_weights <- dbinom(current_obs, size = floor(states$I), prob = fixed_inputs$q, log = TRUE)
  LSE_state_weigths <- matrixStats::logSumExp(log_state_weights)
  param_weight <- exp(LSE_state_weigths)
  if(is.na(param_weight) | param_weight == 0){
    state_norm_weights <- rep(1/length(state_weights), length(state_weights))
  } else {
    state_norm_weights <- exp(log_state_weights - LSE_state_weigths)
  }
  list(state_norm_weights = state_norm_weights, param_weight = param_weight)
}


#' Resample particle indices
#'
#' This function resamples particle \emph{indices} using sampling with replacement.
#' @param n_part Total number of particles, a positive integer.
#' @param norm_weights Normalized weights for all particles, a vector of length
#' \code{n_part}.
#' @return A vector of length \code{n_part}; the element in \code{i}-th position
#' corresponds to the number of times the particle index \code{i} was resampled.
#' @export
branching_sample <- function(n_part, norm_weights){
  sample_res <- sample.int(n = n_part, replace = TRUE, prob = norm_weights)
  as.numeric(table(factor(sample_res, levels = 1:n_part)))
}
