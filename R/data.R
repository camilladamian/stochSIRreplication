#' 7-day rolling mean of Covid-19 daily new cases in Austria
#'
#' A data set created from applying a seven-day rolling mean to Covid-19 daily new cases
#' recorded in Austria from May 1, 2020, to June 15, 2022.
#'
#' @format
#' A data frame with 776 rows and 2 columns:
#' \describe{
#'  \item{date}{Date in YYYY-MM-DD format, an object of class \code{Date}.}
#'  \item{rollmean}{Seven-day rolling mean of Covid-19 daily new cases recorded in Austria.}
#' }
#' @source The original data (daily new cases) were shared under CC-BY-4.0 (\url{https://creativecommons.org/licenses/by/4.0/}) from BMSGPK.
"covid_AT"
