#' effclust (Compute the effective number of clusters)
#'
#' Computes the effective number of clusters in the dataset.
#'
#' @param mx The model object (lm or glm).
#' @param cluster The clustering/grouping variable of interest.
#'
#' @return Outputs a vector of the effective number of clusters per predictor.
#'
#' @examples
#' m1 <- lm(mpg ~ wt, mtcars)
#' effclust(m1, 'cyl')
#'
#' @export
#' @importFrom stats model.matrix nobs
#'
effclust <- function(mx, cluster){

  X <- model.matrix(mx) #design matrix
  n <- nobs(mx) #how many observations
  br <- solve(t(X) %*% X) #bread matrix
  Xj <- split(data.frame(X), cluster) #split by cluster
  js <- lapply(Xj, NROW) #number of clusters per group
  Im <- Yg <- list() #container
  G <- length(Xj) #how many clusters
  k <- ncol(Xj[[1]]) #how many predictors

  # building the i matrix (of 1s)
  for (i in 1:length(js)){
    tmp <- matrix(1, nrow = 1, ncol = js[[i]])
    Im[[i]] <- t(tmp) %*% (tmp)
  }

  ecs <- numeric() #container for effective cluster sizes

  initsel <- rep(0, k) #initialize selection matrix

  for (xx in 1:k){

    sel <- initsel
    sel[xx] <- 1 #which variable to select
    for (i in 1:G){
      Yg[[i]] <- t(sel) %*% br %*% as.matrix(t(Xj[[i]])) %*%
        Im[[i]] %*% as.matrix(Xj[[i]]) %*%
        br %*% sel
    }

    ybar <- mean(unlist(Yg))
    lam <- sum(sapply(Yg, \(x) (x - ybar) / ybar)^2) / G
    ecs[xx] <- G / (1 + lam) #effective cluster size

  }
  names(ecs) <- colnames(X)
  return(ecs)
}
