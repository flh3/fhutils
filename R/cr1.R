#' cr1
#'
#' Returns CR1 standard errors.
#'
#' @param x An lm or glm object.
#' @param dat Original dataset.
#' @param cluster The clustering variable (in quotes).
#' @param ... Extra arguments.
#'
#' @return Returns cluster robust standard errors.
#'
#' @examples
#' # library(sandwich)
#' g1 <- lm(mpg ~ wt + am, data = mtcars)
#' cr1(g1, dat = mtcars, cluster = 'cyl')
#' # sqrt(diag(vcovCL(g1, ~cyl, type = 'HC1')))
#' g2 <- glm(vs ~ wt + am, data = mtcars, family = binomial)
#' cr1(g2, cluster = 'cyl')
#' # sqrt(diag(vcovCL(g2, ~cyl, type = 'HC1')))
#'
#' @importFrom stats residuals vcov
#' @export

cr1 <- function(x, dat = NULL, cluster, ...){

  # if entered without quotes

  # cluster_name <- deparse(substitute(cluster))
  # cluster_name <- gsub('"', '', cluster_name)

  cluster_name = cluster
  ## fast (but not transparent) way of just computing the CRSE
  n <- nobs(x)

  if (is.null(dat)){ #for use in my function
    dat <- eval(x$call$data) #the orig data
  }

  # check if lm or glm
  disp <- if(inherits(x, "glm")) summary(x)$dispersion else summary(x)$sigma^2

  X <- model.matrix(x) #design matrix

  G <- length(unique(dat[[cluster_name]])) #how many clusters
  k <- ncol(X) #how many predictors
  br <- vcov(x) / disp #already computed so fast
  # res <- x$residuals #resid(x)

  res <- if(inherits(x, "glm")) residuals(x, type = "working") * x$weights else
    x$residuals
  u <- X * res #score

  u_g <- rowsum(u, dat[,cluster_name]) #adding them up by group
  df_c <- ((n - 1) / (n - k)) * (G / (G - 1)) #hc1 correction
  M <- (t(u_g) %*% u_g) * df_c #meat matrix

  crse <- sqrt(diag(br %*% M %*% br)) #the CR1 crse
  return(crse)

}




