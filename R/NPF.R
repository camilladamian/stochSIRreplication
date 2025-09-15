#' Obtain posterior estimates through NPF algorithm
#'
#' This function runs the Nested Particle Filter (NPF) algorithm of Crisan and 
#' Míguez (2018) to obtain posterior estimates for states and parameters at each 
#' time point. The function defaults are such that the model of Colaneri et al. 
#' (2025) is assumed.
#'
#' @param n_params Number of parameter particles to sample, a positive scalar. Defaults to 500.
#' @param n_state Number of state particles to sample, a positive scalar. Defaults to 500.
#' @param param_prior_settings Necessary arguments to further specify the chosen 
#' prior distribution for each parameter, a (nested) list. The top-level names must 
#' be \code{k_settings}, \code{s_settings} and \code{m_settings}. In the paper, 
#' the prior distribution of each parameter is chosen to be uniform, so each 
#' top-level element must be a list containing \code{min} and \code{max}. Defaults 
#' to \code{stochSIR_param_prior_settings}, which is available to users upon 
#' package loading.
#' @param state_prior_settings Necessary arguments to further specify the chosen 
#' prior distribution for each state variable, a (nested) list. The top-level 
#' names must be \code{psi_settings} and \code{I_settings}. In the paper, the prior 
#' distribution of \eqn{\psi} is chosen to be normal and the one of \eqn{I} to be 
#' gamma, so the corresponding top-level elements must be a list containing \code{mean} 
#' and \code{sd} and a list containing \code{rate} and \code{shape}. Note that the 
#' paper assumes no-one has recovered at start time, so \eqn{R} is initialized at 
#' 0 for all particles and requires no further specification settings. Defaults to 
#' stochSIR_state_prior_settings, which is available to users upon package loading.
#' @param jitter_settings Necessary arguments to further specify the chosen 
#' jittering kernel, a (nested) list. The top-level names must be \code{k_settings}, 
#' \code{s_settings} and \code{m_settings}. In the paper, the jittering kernel is 
#' chosen to be a truncated Gaussian one and this function calls 
#' \code{truncnorm::rtruncnorm}, so each top-level element must be a list containing 
#' \code{a}, \code{b} and \code{max}. Defaults to \code{stochSIR_jitter_settings}, 
#' which is available to users upon package loading.
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
#' @param tot_time_steps Total number of time steps at which to obtain posterior 
#' estimates, a positive integer. Defaults to \code{length{observations}}.
#' @param delta Time-step size, a positive scalar. Defaults to 1.
#' @param observations (Univariate) time series of observations, a numeric vector. 
#' Defaults to \code{covid_AT$rollmean} (data set \code{covid_AT} is available 
#' to users upon package loading).
#' @export
stochSIR_NPF <- function(n_params = 500, n_state = 500,
                         param_prior_settings = stochSIR_param_prior_settings, 
                         state_prior_settings = stochSIR_state_prior_settings, 
                         jitter_settings = stochSIR_jitter_settings, 
                         fixed_inputs = stochSIR_fixed_inputs, 
                         tot_time_steps = length(observations), delta = 1,
                         observations = covid_AT$rollmean){
  
  ### Initialize
  # Sample from prior
  param_particles <- stochSIR_param_prior(n_params, param_prior_settings)
  state_particles <- vector(mode = "list", length = n_params)
  state_particles <- lapply(state_particles, function(x) stochSIR_state_prior(n_state, state_prior_settings))
  # Pre-allocate outputs list
  outputs <- vector(mode = "list", length = tot_time_steps)
  outputs[[1]] <- stochSIR_NPF_output(state_particles, param_particles, fixed_inputs)
  
  for(i in 2:tot_time_steps){
    previous_obs <- observations[i - 1]
    current_obs <- observations[i]
    
    ### Recursive step
    # Jitter parameter particles
    param_particles <- stochSIR_jitter(n_params, param_particles, jitter_settings)
    # Propagate state particles
    state_particles <- lapply(1:n_params,
                              function(k) stochSIR_state_step(states = state_particles[[k]],
                                                              params = param_particles[k, ],
                                                              fixed_inputs = fixed_inputs, delta = delta,
                                                              previous_obs = previous_obs))
    
    ### Resampling
    # Compute normalized likelihood weights
    likelihood_weights <- lapply(state_particles, 
                                 stochSIR_likelihood_weights, 
                                 current_obs = current_obs, fixed_inputs = fixed_inputs)
    state_norm_weights <- purrr::map(likelihood_weights, "state_norm_weights")
    param_weights <- unlist(purrr::map(likelihood_weights, "param_weight"))
    param_norm_weights <- param_weights/sum(param_weights)
    # Resample with replacement
    state_particles <- lapply(1:n_params, 
                              function(k) state_particles[[k]][rep(1:n_state, times = branching_sample(n_part = n_state, norm_weights = state_norm_weights[[k]])), ])
    params_idx <- rep(1:n_params, 
                      times = branching_sample(n_part = n_params, norm_weights = param_norm_weights))
    param_particles <- param_particles[params_idx, ]
    state_particles <- state_particles[params_idx]
    
    ### Output
    # Compute and format output corresponding to current iteration
    outputs[[i]] <- stochSIR_NPF_output(state_particles, param_particles, fixed_inputs)
  }
  outputs
}