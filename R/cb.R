#' cb
#'
#' Fits a cluster bootstrapped linear model.
#'
#' @param fml A regression formula.
#' @param cluster The clustering variable (in quotes).
#' @param dat Original dataset.
#' @param B Number of bootstrap replicates.
#' @param group Specify a group (in quotes) to stratify on.
#' @param interval Type of confidence interval ('bca' = bias corrected and accelerated or 'percentile').
#' @param seed A seed for replicability.
#'
#' @return Regression results.
#'
#' @examples
#' g1 <- cb(mpg ~ wt + am, cluster = 'cyl', dat = mtcars)
#
#'
#' @importFrom stats coef lm model.frame quantile update as.formula lm.fit qnorm pnorm sd
#' @export
cb <- function(fml, cluster, dat, B = 1000, group = NULL, interval = "bca", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Initial Setup
  mf <- model.frame(fml, data = dat)
  X_orig <- model.matrix(fml, mf)
  y_orig <- model.response(mf)

  cluster_var <- as.character(dat[[cluster]])
  cluster_ids <- unique(cluster_var)
  G <- length(cluster_ids)

  # 2. Optional Grouping (Stratification)
  if (!is.null(group)) {
    group_var   <- as.character(dat[[group]])
    cluster_group_counts <- tapply(group_var, cluster_var, function(x) length(unique(x)))
    if (any(cluster_group_counts > 1)) {
      stop("Error: Some clusters are associated with multiple groups.")
    }
    cluster_map <- unique(data.frame(cluster = cluster_var, group = group_var))
    clusters_by_group <- split(cluster_map$cluster, cluster_map$group)
  } else {
    clusters_by_group <- list(all = cluster_ids)
  }

  # Optimization: Pre-split data
  X_list <- lapply(split(as.data.frame(X_orig), cluster_var), as.matrix)
  y_list <- split(y_orig, cluster_var)

  k <- ncol(X_orig)
  boot_coefs <- matrix(NA, nrow = B, ncol = k)
  colnames(boot_coefs) <- colnames(X_orig)

  # 3. Bootstrap Loop
  for (b in 1:B) {
    resampled_clusters <- unlist(lapply(clusters_by_group, function(cls) {
      sample(cls, size = length(cls), replace = TRUE)
    }), use.names = FALSE)

    X_star <- do.call(rbind, X_list[resampled_clusters])
    y_star <- unlist(y_list[resampled_clusters], use.names = FALSE)
    boot_coefs[b, ] <- lm.fit(X_star, y_star)$coefficients
  }

  # 4. Original Fit
  # orig_fit <- lm.fit(X_orig, y_orig)
  # orig_coefs <- orig_fit$coefficients

  orig_fit <- lm(fml, data = dat)
  orig_coefs <- coef(orig_fit)
  orig_se <- sqrt(diag(vcov(orig_fit)))


  # 5. Confidence Interval Logic
  ci_low <- numeric(k)
  ci_high <- numeric(k)

  if (interval == "percentile") {
    ci_bounds <- apply(boot_coefs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    ci_low <- ci_bounds[1, ]
    ci_high <- ci_bounds[2, ]
  } else if (interval == "bca") {
    message("Calculating BCa Acceleration via Cluster Jackknife...")

    # Cluster Jackknife for Acceleration (a)
    jk_coefs <- matrix(NA, nrow = G, ncol = k)
    for (i in 1:G) {
      # Identify all indices EXCEPT those in the current cluster
      keep_idx <- which(cluster_var != cluster_ids[i])
      jk_fit <- lm.fit(X_orig[keep_idx, , drop = FALSE], y_orig[keep_idx])
      jk_coefs[i, ] <- jk_fit$coefficients
    }

    # Calculate BCa for each parameter
    for (j in 1:k) {
      # Bias correction (z0)
      # z0 = qnorm(proportion of boot estimates < original estimate)
      p_less <- mean(boot_coefs[, j] < orig_coefs[j], na.rm = TRUE)
      z0 <- qnorm(p_less)

      # Acceleration (a) using Jackknife estimates
      jk_mean <- mean(jk_coefs[, j])
      diffs <- jk_mean - jk_coefs[, j]
      accel <- sum(diffs^3) / (6 * (sum(diffs^2))^1.5)

      # Adjusted quantiles
      z_alpha <- qnorm(c(0.025, 0.975))
      adj_probs <- pnorm(z0 + (z0 + z_alpha) / (1 - accel * (z0 + z_alpha)))

      # Final CI from bootstrap distribution
      ci_val <- quantile(boot_coefs[, j], probs = adj_probs, na.rm = TRUE)
      ci_low[j] <- ci_val[1]
      ci_high[j] <- ci_val[2]
    }
  }

  # 6. Final Results
  boot_se <- apply(boot_coefs, 2, sd, na.rm = TRUE)
  res <- data.frame(
    term = names(orig_coefs),
    estimate = as.numeric(orig_coefs),
    orig.std.error = round(as.numeric(orig_se), 8),
    boot.std.error = as.numeric(boot_se),
    statistic = as.numeric(orig_coefs / boot_se),
    conf.low = ci_low,
    conf.high = ci_high
  )

  cat("\nResampling units:", cluster)
  if(!is.null(group)) cat(" | Stratified by:", group)
  cat(" | Interval Type:", interval)
  cat(" | Iterations:", B, "\n\n")

  return(res)
}
