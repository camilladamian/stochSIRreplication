#' Jitter parameter particles
#'
#' This function 'jitters' \emph{all} parameter particles \eqn{\theta = (\kappa, \sigma, \mu)} 
#' using a truncated Gaussian kernel with mean corresponding to the currently 
#' available parameter particles and fixed variance (see Crisan and Míguez (2018) 
#' for more details about the possible choices of jittering kernels).
#'
#' @param n_params Number of parameter particles to jitter, a positive scalar.
#' @param param_particles Parameter particles to jitter, a data frame 
#' with columns \code{k}, \code{s} and \code{m}.
#' @param jitter_settings Necessary arguments to further specify the chosen 
#' jittering kernel, a (nested) list. The top-level names must be \code{k_settings}, 
#' \code{s_settings} and \code{m_settings}. In the paper, the jittering kernel is 
#' chosen to be a truncated Gaussian one and this function calls 
#' \code{truncnorm::rtruncnorm}, so each top-level element must be a list containing 
#' \code{a}, \code{b} and \code{max}. Defaults to \code{stochSIR_jitter_settings}, 
#' which is available to users upon package loading.
#' @return A data frame with \code{n_params} rows and columns \code{k}, \code{s} 
#' and \code{m} (corresponding to \eqn{\kappa, \sigma} and \eqn{\mu}, respectively).
#' @export
stochSIR_jitter <- function(n_params, param_particles, 
                            jitter_settings = stochSIR_jitter_settings){
  data.frame(k = do.call(truncnorm::rtruncnorm, c(n = n_params, list(mean = param_particles$k), 
                                                  jitter_settings$k_settings)),
             s = do.call(truncnorm::rtruncnorm, c(n = n_params, list(mean = param_particles$s), 
                                                  jitter_settings$s_settings)),
             m = do.call(truncnorm::rtruncnorm, c(n = n_params, list(mean = param_particles$m), 
                                                  jitter_settings$m_settings)))
}

#' Propagate state particles from their transition kernel
#'
#' For a given parameter particle, this function propagates state particles \eqn{x = (\psi, i, r)}
#' from their transition kernel.
#'
#' @param states A data frame with columns \code{psi}, \code{I} and \code{R} 
#' (corresponding to \eqn{\psi, i} and \eqn{r}, respectively).
#' @param params A named vector with elements \code{k}, \code{m} and \code{s} 
#' (corresponding to \eqn{\kappa, \sigma} and \eqn{\mu}, respectively).
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
#' @param delta Time-step size, a positive scalar.
#' @param previous_obs Number of people detected at the previous time point and 
#' removed from pool of infectious people, a positive integer. Only used if 
#' \code{quarantine} is \code{TRUE}.
#' @return A data frame with the same number of rows as argument \code{states} and 
#' columns \code{psi}, \code{I} and \code{R} (corresponding to \eqn{\psi, i} and 
#' \eqn{r}, respectively).
#' @export
stochSIR_state_step <- function(states, params, fixed_inputs = stochSIR_fixed_inputs, delta, previous_obs){
  list2env(as.list(params), envir = environment())
  list2env(as.list(states), envir = environment())
  list2env(as.list(fixed_inputs), envir = environment())
  psi <- psi + k * (m - psi) * delta + s * sqrt(delta) * rnorm(length(psi))
  S <- tot_pop - (I + R)
  intensity <- I * exp(psi) * (S / tot_pop)
  I_plus <- rpois(length(intensity), intensity * delta)
  if(quarantine){
    I_minus <- previous_obs + g * I
  } else {
    I_minus <- g * I
  }
  I <- I + I_plus - I_minus
  R <- (1 - d_f) * R + I_minus
  I[I < 0] <- 0
  data.frame(psi = psi, I = I, R = R)
}