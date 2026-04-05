wcr <- function(fmlf, dat, param, B = 999, cluster){

  # fmlr <- 'wage ~ ttl_exp + collgrad'
  # fmlf <- 'wage ~ tenure + ttl_exp + collgrad'

  fmlr <- update(fmlf, as.formula(paste("~ . -", param)))
  # 1. Setup and Restricted Model (WCR)
  # Testing the hypothesis that the coefficient of 'wt' is 0

  orig <- lm(fmlf, data = dat)
  gr <- dat[[cluster]]

  X <- model.matrix(orig) # Full X matrix [cite: 188]
  se <- cr1(orig, dat = dat, cluster = cluster)
  tobs <- coef(orig) / se
  G <- length(unique(gr))
  B <- 9999
  k <- ncol(X)
  n <- nrow(X)
  # gr <- dat[[cluster]]


  g_res <- lm(fmlr, data = dat) # Null imposed: 'wt' removed [cite: 193]

   utilde <- resid(g_res) # get residuals of restricted


  # 2. Pre-calculate Constant "Bread" and Small Sample Adjustment [cite: 127, 524, 646]
  invXX <- solve(crossprod(X))
  adj <- (G / (G - 1)) * ((n - 1) / (n - k))

  # 3. Predict under the Null to get the center of the bootstrap DGP [cite: 163, 169]
  pred_null <- predict(g_res, newdata = dat)

  # 4. Collapse X and Residuals into Cluster-Level "Scores" [cite: 548, 610]
  # This is the O(N) step performed only once [cite: 614]
  SX_u <- rowsum(X * utilde, gr) # G x k matrix of cluster-scores


  # S <- Matrix::fac2sparse(gr) #same as above
  # scores <- X * as.vector(utilde)
  # U_g <- S %*% scores #faster

  # 5. Generate Wild Weights (Rademacher) [cite: 166, 274, 649]
  v <- matrix(sample(c(-1, 1), G * B, replace = TRUE), nrow = G, ncol = B)

  # 6. FAST NUMERATOR: R(beta_star - beta_null) [cite: 564, 655]
  # Directly compute the k x B matrix of coefficient deviations
  # (invXX %*% t(X)) %*% (utilde * v[at cluster level]) simplifies to:
  numer_all <- invXX %*% t(SX_u) %*% v # k x B matrix [cite: 564]

  # 7. FAST DENOMINATOR: The "Meat" for each B [cite: 607, 610]
  # We only need the variance of the specific coefficient (e.g., 'wt' is index 2)
  # The paper uses a G x G projection to avoid N x N matrices [cite: 610, 613]
  # For R, the most memory-efficient way is to loop over the B iterations
  # of the G x k scores:
  weights_v <- v  # G x B


  boot_se <- vapply(1:B, function(b) {
    # Perturb the scores directly by the wild weights
    scores_b <- SX_u * weights_v[, b] # G x k
    meat_b <- crossprod(scores_b)    # k x k
    vcov_b <- invXX %*% meat_b %*% invXX * adj
    sqrt(diag(vcov_b))
  }, numeric(k)) # k x B matrix of SEs [cite: 196]



  # 8. P-Value Calculation [cite: 203, 472, 662]
  t_boot <- numer_all / boot_se # Resulting t-stats for all B
  # t_orig <- coef(summary(lm_robust(as.formula(fmlf), se_type = 'stata', data = data, cluster = industry)))["tenure", "t value"]
  t_orig <- tobs[param]

  p_val <- mean(abs(t_boot[param, ]) >= abs(t_orig)) # Symmetric p-value [cite: 203, 821]
  cat("p-value:", p_val, "\n")
}
