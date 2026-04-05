#' cv (Coefficient of Variation)
#'
#' Computes the coefficient of variation.
#'
#' @param cluster The clustering/grouping variable of interest.
#'
#' @return Outputs the coefficient of variation which represents the degree of
#' cluster size imbalance (0 = completely balanced).
#'
#' @examples
#' cv(mtcars$cyl)
#'
#' @export
#'
cv <- function(cluster){
  xx <- table(cluster)
  ##nm <- names(xx)
  cv <-  sd(xx) / mean(xx)
  return(cv)
}
