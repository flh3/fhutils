wcr_optimized <- function(fmlf, dat, param, B = 999, cluster) {

  # 1. Setup and Restricted Model (WCR)
  fmlr <- update(fmlf, as.formula(paste("~ . -", param)))

  orig <- lm(fmlf, data = dat)
  X <- model.matrix(orig)
  k <- ncol(X)
  N <- nrow(X)

  gr <- dat[[cluster]]
  G <- length(unique(gr))

  # Small-sample adjustment factor (m)
  m <- (G / (G - 1)) * ((N - 1) / (N - k))

  g_res <- lm(fmlr, data = dat)
  utilde <- as.numeric(resid(g_res))

  # 2. Extract the parameter restriction into a conformable column vector R0
  param_idx <- which(colnames(X) == param)
  if (length(param_idx) == 0) stop("Parameter not found in the model.")
  R0_col <- matrix(0, nrow = k, ncol = 1)
  R0_col[param_idx, 1] <- 1

  # 3. Pre-calculate Constant "Bread"
  invXX <- solve(crossprod(X))

  # 4. Generate Rademacher Weights for all replications at once
  # We use G x (B+1) and set the first column to 1s to perfectly reproduce
  # the original sample's t-statistic alongside the bootstrap stats.
  v <- matrix(sample(c(-1, 1), G * (B + 1), replace = TRUE), nrow = G)
  v[, 1] <- 1

  # ---------------------------------------------------------------------
  # FAST EXECUTION MATRICES: Dimension reduction from N to G
  # ---------------------------------------------------------------------

  # Multiply thin matrix first: X * (X'X)^-1 * R'
  XinvXXR <- X %*% (invXX %*% R0_col)

  # FAST NUMERATOR COMPONENTS: Collapse to G x 1 cluster-scores before the B loop
  SXinvXXRu <- rowsum(XinvXXR * utilde, gr)

  # FAST DENOMINATOR COMPONENTS: Collapse X and residuals to cluster-scores
  SX_XinvXXR <- rowsum(X * as.vector(XinvXXR), gr)
  SX_utilde <- rowsum(X * utilde, gr)

  # The fixed 'J' Matrix component (G x G)
  J_matrix <- diag(as.vector(SXinvXXRu)) - (SX_XinvXXR %*% invXX %*% t(SX_utilde))

  # ---------------------------------------------------------------------
  # VECTORIZED COMPUTATION FOR ALL REPLICATIONS SIMULTANEOUSLY
  # ---------------------------------------------------------------------

  # Compute k x B Numerators (O(GB) operation)
  numer <- t(SXinvXXRu) %*% v

  # Compute B Denominators (O(G^2 B) operation)
  J <- J_matrix %*% v

  # Resulting t-statistics for actual sample (position 1) and all B replications
  t_stats <- abs(numer) / t(sqrt(m * colSums(J^2)))

  # 8. P-Value Calculation
  t_orig <- t_stats[9]       # The actual sample's cluster-robust t-statistic
  t_boot <- t_stats[-1]      # The B bootstrap t-statistics

  p_val <- sum(t_boot > t_orig) / B

  cat("p-value:", p_val, "\n")
  return(p_val)
}
