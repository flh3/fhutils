#' ivex (instrumental variable example)
#'
#' Simulated data to demonstrate instrumental variable estimation.
#'
#' @format A data frame with 500 rows and 5 variables:
#' \describe{
#'   \item{u}{A confounder.}
#'   \item{z}{Treatment assignment dummy code (1 = treatment, 0 = control).}
#'   \item{tr}{Treatment takeup dummy code (1 = took the treatment; 0 = did not take treatment).}
#'   \item{y}{The outcome of interest.}
#'   \item{x1}{A covariate.}
#' }
"ivex"
