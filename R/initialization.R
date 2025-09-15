#' Sample parameter particles from their prior distribution
#'
#' This function samples parameter particles \eqn{\theta = (\kappa, \sigma, \mu)}
#' from their prior distribution.
#'
#' @param n_params Number of parameter particles to sample, a positive scalar.
#' @param param_prior_settings Necessary arguments to further specify the chosen 
#' prior distribution for each parameter, a (nested) list. The top-level names must 
#' be \code{k_settings}, \code{s_settings} and \code{m_settings}. In the paper, 
#' the prior distribution of each parameter is chosen to be uniform, so each 
#' top-level element must be a list containing \code{min} and \code{max}. Defaults 
#' to \code{stochSIR_param_prior_settings}, which is available to users upon 
#' package loading.
#' @return A data frame with \code{n_params} rows and columns \code{k}, \code{s} 
#' and \code{m} (corresponding to \eqn{\kappa, \sigma} and \eqn{\mu}, respectively).
#' @export
stochSIR_param_prior <- function(n_params, 
                                 param_prior_settings = stochSIR_param_prior_settings){
  data.frame(k = do.call(runif, c(n = n_params, param_prior_settings$k_settings)),
             s = do.call(runif, c(n = n_params, param_prior_settings$s_settings)),
             m = do.call(runif, c(n = n_params, param_prior_settings$m_settings))
  )
}

#' Sample state particles from their prior distribution
#'
#' This function samples state particles \eqn{x = (\psi, i, r)}
#' from their prior distribution.
#'
#' @param n_state Number of state particles to sample, a positive scalar.
#' @param state_prior_settings Necessary arguments to further specify the chosen 
#' prior distribution for each state variable, a (nested) list. The top-level 
#' names must be \code{psi_settings} and \code{I_settings}. In the paper, the prior 
#' distribution of \eqn{\psi} is chosen to be normal and the one of \eqn{I} to be 
#' gamma, so the corresponding top-level elements must be a list containing \code{mean} 
#' and \code{sd} and a list containing \code{rate} and \code{shape}. Note that the 
#' paper assumes no-one has recovered at start time, so \eqn{R} is initialized at 
#' 0 for all particles and requires no further specification settings. Defaults to 
#' stochSIR_state_prior_settings, which is available to users upon package loading.
#' @return A data frame with \code{n_state} rows and columns \code{psi}, \code{I} 
#' and \code{R} (corresponding to \eqn{\psi, i} and \eqn{r}, respectively).
#' @export
stochSIR_state_prior <- function(n_state, 
                                 state_prior_settings = stochSIR_state_prior_settings){
  data.frame(psi = do.call(rnorm, c(n = n_state, state_prior_settings$psi_settings)),
             I = floor(do.call(rgamma, c(n = n_state, state_prior_settings$I_settings))),
             R = rep(0, n_state))
}