#' r2mm (Computes r2 for multilevel models)
#'
#' Computes the r2 for a multilevel model. Provides marginal and conditional r2.
#'
#' @param obj Output using lmer (merMod)
#'
#' @return A data frame
#'
#' @examples
#' library(lme4)
#' data(pkrdd)
#' m1 <- lmer(abcs ~ avar * takeup + female + (1 | tr), data = pkrdd)
#' r2mm(m1)
#'
#' @importFrom lme4 fixef VarCorr ranef getME
#' @importFrom stats cov sigma
#' @importFrom Matrix bdiag
#'
#' @export
#'
r2mm <- function(obj){

  X <- model.matrix(obj)
  b <- fixef(obj)
  sig2 <- sigma(obj)^2
  vc2 <- VarCorr(obj)
  vc <- Matrix::bdiag(vc2)
  cvf <- cov(X)
  fv <- as.numeric(t(b) %*% cvf %*% b)
  clust <- names(ranef(obj))

  rf0 <- ranef(obj)
  levs <- length(rf0)
  nore <- sum(sapply(rf0, length)) #number of RE

  #cat("Number of random effects:", nore)
  cat('\nNumber of levels:', levs + 1, '\n')
  nors <- nore - levs #number of random slopes

  if (levs == nore){
    cat('Type of Model: Random Intercept Model \n')
    tau00 = sum(diag(vc))
    rsv = 0

  } else {

    cat('Type of Model: Random Slope Model \n')
    rsdat <- getME(obj, 'mmList')
    tx <- list()
    for (i in 1:levs){
      m <- colMeans(rsdat[[i]])
      tx[[i]] <- as.numeric(t(m) %*% vc2[[i]] %*% m)
    }

    tau00 <- sum(unlist(tx)) #random intercept variance
    # print(tau00)

    trs <- list()
    for (i in 1:levs){
      rsl <- ncol(rsdat[[i]])

      if (rsl == 1){
        trs[[i]] <- 0
      } else {
        vct <- bdiag(vc2[[i]])[2:rsl, 2:rsl]
        cv2 <- cov(rsdat[[i]][, -1, drop = F])
        trs[[i]] <- sum(diag(cv2 %*% vct))
      }
    }

    rsv <- sum(unlist(trs))
  }


  totv <- fv + tau00 + rsv + sig2
  fv2 <- fv / totv
  rs2 <- rsv / totv
  ri2 <- tau00 / totv
  si2 <- sig2 / totv
  con <- (fv2 + rs2 + ri2)

  data.frame(marginal = fv2, conditional = con, rs = rs2, ri = ri2, sigma2 = si2,
             totv = totv)
}
