<div align="justify">

# Introduction
The replication package for the empirical analysis included in our paper *A Filtering 
Approach for Statistical Inference in a Stochastic SIR Model with an Application 
to Covid-19 Data* can be installed (if needed) from GitHub and loaded as follows:
```r
### Install if needed (do this exactly once)
# devtools::install_github("camilladamian/stochSIRreplication")
### Load
library(stochSIRreplication)
```
### Note
The package `stochSIRreplication` only includes functions and data meant to replicate our 
empirical analysis, and it is not a package aimed at a general implementation of 
stochastic SIR models and at their estimation. In particular, the transition and 
observation densities are assumed to be those of Section 2 of our paper, and the 
default values for fixed parameters are set according to those reported in Appendix 
C of the Supplementary Material accompanying our paper.

# Data Description and Visualization
In the paper *A Filtering Approach for Statistical Inference in a Stochastic SIR 
Model with an Application to Covid-19 Data*, we apply the Nested Particle Filter 
(NPF) approach of Crisan and Míguez (2018)[^1] to Austrian Covid-19 data from May 1, 
2020 to June 15, 2022. In particular, we use a seven-day rolling average of confirmed 
cases to avoid weekly seasonality effects, such as the fewer tests performed over 
the weekend. The data set `covid_AT` is available to users upon loading the package, 
and the data can be visualized as follows:

```r
plot(covid_AT$date, covid_AT$rollmean, type = "l", col = "slateblue", xlab = "",
     ylab = "", lwd = 2, main = "Austria: Seven-day Rolling Average of
     Positive Tests")
```
<img width="960" height="576" alt="ATdata" src="https://github.com/user-attachments/assets/44b0c0bd-f005-41a4-a147-1a0284998946" />

[^1]: Crisan, D. and Miguez, J. (2018). Nested particle filters for online parameter estimation in 
discrete-time state space Markov models. *Bernoulli* **24**, 3039–3086.

# Nested Particle Filter (NPF)
To replicate the results of our paper *A Filtering Approach for Statistical Inference 
in a Stochastic SIR Model with an Application to Covid-19 Data*, one can simply 
use the `stochSIR_NPF` function with all arguments set to their defaults as follows 
(please note that this might take around 6.5 minutes to run):
```r
set.seed(2727)
stochSIR_NPF_results <- stochSIR_NPF()
```

For a detailed explanation of the function arguments and their defaults, please 
refer to Appendix C of the Supplementary Material accompanying our paper and/or 
to the function documentation (available by running `?"stochSIR_NPF"`).

The object `stochSIR_NPF_results` obtained above is then a list-of-lists of length 
776. To extract all relevant quantities in a usable format, please use:
```r
rep_rate_res <- data.frame(do.call(rbind, purrr::map(stochSIR_NPF_results , "rep_rate_res")))
k_res <- data.frame(do.call(rbind, purrr::map(stochSIR_NPF_results , "k_res")))
s_res <- data.frame(do.call(rbind, purrr::map(stochSIR_NPF_results , "s_res")))
m_res <- data.frame(do.call(rbind, purrr::map(stochSIR_NPF_results , "m_res")))
ratio_res <- data.frame(do.call(rbind, purrr::map(stochSIR_NPF_results , "ratio_res")))
```
# NPF Results Visualization
### Effective Reproduction Rate
```{r, echo = TRUE, fig.width = 10, fig.height = 6, out.width = "75%", out.height = "75%", fig.align = "center"}
plot(covid_AT$date, rep_rate_res$rep_rate_hat, type = "l", xlab = "", ylab = "", 
     main = "Filtered Effective Reproduction Rate", lwd = 2, col = "coral")
```
<img width="960" height="576" alt="rep_rate" src="https://github.com/user-attachments/assets/9079cc11-a4e4-454b-bf60-24524c22cdaf" />

### Parameter Posterior Evolution Over Time
```{r, fig.width = 15, fig.height = 4, out.width = "75%", out.height = "75%", fig.align = "center"}
op <- par(mfrow = c(1, 3))
# Mu: posterior evolution
with(m_res, plot(covid_AT$date, m_lq, type = "l", ylim = range(c(m_lq, m_uq)), 
                 col = "slateblue", lwd = 2, xlab = "", ylab = "",
                 main = expression(paste(mu, ": Parameter Estimation")), 
                 cex.main = 2, cex.axis = 1.5, cex.lab = 1.5))
with(m_res, lines(covid_AT$date, m_uq, col = "slateblue", lwd = 2))
with(m_res, lines(covid_AT$date, m_hat, col = "coral", lwd = 2))
legend("topright", legend = c("95% quantile", "Mean", "5% quantile"),  bty = "n", 
       col = c("slateblue", "coral", "slateblue"), lwd = rep(2, 4), cex = 1.25)

# Kappa: posterior evolution
with(k_res, plot(covid_AT$date, k_lq, type = "l", ylim = range(c(k_lq, k_uq)), 
                 col = "slateblue", lwd = 2, xlab = "", ylab = "",
                 main = expression(paste(kappa, ": Parameter Estimation")), 
                 cex.main = 2, cex.axis = 1.5, cex.lab = 1.5))
with(k_res, lines(covid_AT$date, k_uq, col = "slateblue", lwd = 2))
with(k_res, lines(covid_AT$date, k_hat, col = "coral", lwd = 2))
legend("topright", legend = c("95% quantile", "Mean", "5% quantile"),  bty = "n", 
       col = c("slateblue", "coral", "slateblue"), lwd = rep(2, 4), cex = 1.25)

# Sigma: posterior evolution
with(s_res, plot(covid_AT$date, s_lq, type = "l", ylim = range(c(s_lq, s_uq)), 
                 col = "slateblue", lwd = 2, xlab = "", ylab = "",
                 main = expression(paste(sigma, ": Parameter Estimation")), 
                 cex.main = 2, cex.axis = 1.5, cex.lab = 1.5))
with(s_res, lines(covid_AT$date, s_uq, col = "slateblue", lwd = 2))
with(s_res, lines(covid_AT$date, s_hat, col = "coral", lwd = 2))
legend("topright", legend = c("95% quantile", "Mean", "5% quantile"),  bty = "n", 
       col = c("slateblue", "coral", "slateblue"), lwd = rep(2, 4), cex = 1.25)
par(op)
```
<img width="1440" height="384" alt="posteriors" src="https://github.com/user-attachments/assets/f350f614-915c-44b8-a290-a15cb0164a5d" />

</div>
