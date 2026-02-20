#' pkrdd (Preschool data for regression discontinuity designs)
#'
#' Simulated data to test out regression discontinuity designs and instrumental variables.
#'
#' @format A data frame with 161 rows and 8 variables:
#' \describe{
#'   \item{abcs}{Outcome of interest. The number of letters known.}
#'   \item{avar}{The assignment (or forcing) variable (based on date of birth).}
#'   \item{takeup}{If the student followed the treatment assignment.}
#'   \item{tr}{If the child should be in kindergarten (1) or preschool (0)}
#'   \item{white}{If the student race was White (1) or nonWhite(0).}
#'   \item{dis}{Disability status (1 = yes / 0 = no).}
#'   \item{female}{Whether student was female (1) or male (0).}
#'   \item{sch2}{The school the child attended.}
#' }
#' @source \url{https://www.sciencedirect.com/science/article/abs/pii/S088520061630062X}
"pkrdd"
