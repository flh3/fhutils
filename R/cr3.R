#' cr3
#'
#' Returns CR3 (jackknifed) standard errors.
#'
#' @param model An lm object.
#' @param dat Original dataset.
#' @param cluster The clustering variable (in quotes).
#' @param df Specify degrees of freedom.
#' @param adj Includes few cluster adjustments.
#' @param seonly Show the standard error only.
#'
#' @return Returns cluster robust standard errors (seonly = TRUE) or
#' complete regression output.
#'
#' @examples
#' g1 <- lm(mpg ~ wt + am, data = mtcars)
#' cr3(g1, dat = mtcars, cluster = 'cyl')
#
#'
#' @importFrom stats residuals vcov coef model.response pt
#' @export
#'
cr3 <- function(model, dat = NULL, cluster, df = Inf, adj = TRUE, seonly = TRUE) {

  # 1. Extract the data actually used in the model (handles NAs automatically)
  mod_data <- model$model
  X <- model.matrix(model)
  y <- model.response(mod_data)
  cf <- coef(model)
  e <- model$residuals
  nm <- names(coef(model))
  # cluster_col <- deparse(substitute(cluster))
  # cluster_col <- gsub('["\']', '', cluster_col) # Clean quotes

  cluster_col = cluster

  if (is.null(dat)) dat <- eval(model$call$data)

  if (!is.null(attr(mod_data, "na.action"))) {
    rows_kept <- setdiff(seq_len(nrow(dat)), attr(mod_data, "na.action"))
    cluster_vec <- as.character(dat[rows_kept, cluster_col])
  } else {
    cluster_vec <- as.character(dat[[cluster_col]])
  }

  # 3. Pre-calculate (X'X)^-1
  xtx_inv <- chol2inv(chol(crossprod(X)))

  # 4. Grouping
  idx_list <- split(seq_len(nrow(X)), cluster_vec)
  n_clusters <- length(idx_list)

  # 5.
  diff_sq <- matrix(0, nrow = n_clusters, ncol = length(cf))

  for (i in seq_along(idx_list)) {
    ii <- idx_list[[i]]
    Xg <- X[ii, , drop = FALSE]
    eg <- e[ii]

    # Calculate the cluster leverage matrix
    Hg <- Xg %*% xtx_inv %*% t(Xg)

    # Identify size of the identity matrix for this cluster
    Ig <- diag(length(ii))

    # Solve (I - Hg)^-1 %*% eg
    # This is the "Leave-one-cluster-out" adjustment
    adj_e <- solve(Ig - Hg, eg)

    # Compute the change in coefficients: (X'X)^-1 %*% Xg' %*% adjusted_errors
    bg_diff <- xtx_inv %*% t(Xg) %*% adj_e
    diff_sq[i, ] <- as.numeric(bg_diff)^2
  }

  # 6. Final Stats
  vc_diag <- colSums(diff_sq)
  adjf <- if(adj) (n_clusters - 1) / n_clusters else 1
  se <- sqrt(adjf * vc_diag)
  names(se) <- nm

  tstat <- cf / se
  if (is.infinite(df)) df <- n_clusters - 1
  pv <- 2 * pt(-abs(tstat), df)

  if (seonly) return(se = se)
  else
    return(data.frame(
      est = cf, se = se, t = tstat, df = df, pv = pv
    ))
}


