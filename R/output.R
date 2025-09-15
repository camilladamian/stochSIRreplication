#' Compute and format relevant NPF outputs
#'
#' This function computes relevant characteristics of the approximate posterior 
#' distribution and formats them 'by variable'.
#'
#' @param state_particles List of state particles, each element of which is a data 
#' frame with columns \code{psi}, \code{I} and \code{R} (corresponding to \eqn{\psi, i} 
#' and \eqn{r}, respectively).
#' @param param_particles Parameter particles, a data frame 
#' with columns \code{k}, \code{s} and \code{m}.
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
#' @return A list with top-level names \code{rep_rate_res}, \code{k_res}, 
#' \code{s_res}, \code{m_res} and \code{ratio_res} (corresponding to the reproduction 
#' number \eqn{\mathcal{R}}, the parameters \eqn{\kappa}, \eqn{\sigma}, \eqn{\mu} 
#' and the ratio \eqn{\sigma^2/2\kappa}, respectively). Each top-level element is 
#' a data frame containing the mean, the 5% and the 95% quantile of the posterior 
#' distribution of each variable of interest.
#' @export
stochSIR_NPF_output <- function(state_particles, param_particles, fixed_inputs){
  list2env(as.list(fixed_inputs), envir = environment())
  psi_particles <- unlist(purrr::map(state_particles, "psi"))
  I_particles <- unlist(purrr::map(state_particles, "I"))
  R_particles <- unlist(purrr::map(state_particles, "R"))
  k_particles <- param_particles$k
  m_particles <- param_particles$m
  s_particles <- param_particles$s
  beta_particles <- exp(psi_particles)
  S_particles <- tot_pop - (I_particles + R_particles)
  rep_rate_particles <- c(beta_particles/(g + q - g*q) * S_particles/tot_pop)
  ratio_particles <- (s_particles^2)/(2 * k_particles)
  rep_rate_res <- c(rep_rate_hat = mean(rep_rate_particles), 
                    rep_rate_lq = quantile(rep_rate_particles, 0.05, names = FALSE), 
                    rep_rate_uq = quantile(rep_rate_particles, 0.95, names = FALSE))
  k_res <- c(k_hat = mean(k_particles), 
             k_lq = quantile(k_particles, 0.05, names = FALSE), 
             k_uq = quantile(k_particles, 0.95, names = FALSE))
  s_res <- c(s_hat = mean(s_particles), 
             s_lq = quantile(s_particles, 0.05, names = FALSE), 
             s_uq = quantile(s_particles, 0.95, names = FALSE))
  m_res <- c(m_hat = mean(m_particles), 
             m_lq = quantile(m_particles, 0.05, names = FALSE), 
             m_uq = quantile(m_particles, 0.95, names = FALSE))
  ratio_res <- c(ratio_hat = mean(ratio_particles), 
                 ratio_lq = quantile(ratio_particles, 0.05, names = FALSE), 
                 ratio_uq = quantile(ratio_particles, 0.95, names = FALSE))
  list(rep_rate_res = rep_rate_res, k_res = k_res, s_res = s_res, 
       m_res = m_res, ratio_res = ratio_res)
}