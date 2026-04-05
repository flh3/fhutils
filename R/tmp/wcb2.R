wcb2 <- function(fml, cluster, dat, B = 999, seed = 0, cr3 = FALSE,
          param = NULL)
{
  if (!is.null(param)) {
    m1 <- lm(fml, data = dat)
    n <- nobs(m1)
    Xo <- model.matrix(m1)
    # obread <- vcov(m1)/(summary(m1)$sigma^2)
    obread <- solve(crossprod(Xo))
    if (!cr3) {
      obsse2 <- cr1(m1, cluster = cluster, dat = dat)


    }
    else {
      obsse2 <- cr3(m1, cluster = cluster, dat = dat)


    }
    tb2 = coef(m1)[param]
    tstat = tb2/obsse2[param]
    drop_var <- param
    fml <- update(fml, as.formula(paste("~ . -", drop_var)))
  }
  m1 <- lm(fml, data = dat)
  tb <- coef(m1)
  gr <- dat[[cluster]]
  n <- nobs(m1)
  y <- model.frame(m1)[[1]]
  yhat <- m1$fitted.values
  re <- m1$residuals
  clust_ids <- unique(gr)
  G <- length(clust_ids)
  k <- length(tb)
  X <- model.matrix(m1)
  # bread <- vcov(m1)/(summary(m1)$sigma^2)
  bread <- solve(crossprod(X))
  # if (!is.null(param))
  #   k = k + 1
  df_c <- ((n - 1)/(n - k)) * (G/(G - 1))


  if (!cr3) { #cr1
    # obsse <- cr1(m1, cluster = cluster, dat = dat)

    u_g <- rowsum(X * as.numeric(re), gr)
    meat <- (t(u_g) %*% u_g) * df_c

    # 3. Assemble the Sandwich: Bread %*% Meat %*% Bread
    vcov <- bread %*% meat %*% bread
    obsse <- sqrt(diag(vcov))
  }
  else { #cr3
    obsse <- cr3(m1, cluster = cluster, dat = dat)
  }
  tt <- tb/obsse

  tstar <- matrix(NA, ncol = k, nrow = B)
  tsr <- numeric()



  #######

  if (seed != 0)
    set.seed(seed)

  W_matrix <- matrix(sample(c(-1, 1), G * B, replace = TRUE),
                     nrow = G, ncol = B)
  W_N <- W_matrix[match(gr, clust_ids), ]
  Y_star <- yhat + (re * W_N)
  B_star <- bread %*% t(X) %*% Y_star

  print(B_star)
  rb <<- Y_star - X %*% B_star
  X <<- X
  print(dim(rb))

  U_g = rowsum(X * rb, gr)
  print(U_g)

  boot_se_matrix <- matrix(NA, nrow = k, ncol = B)

  for (b in 1:B) {
    # 1. Get residuals for this bootstrap iteration
    # resid = y_star - X %*% beta_star
    rb <- Y_star[, b] - (X %*% B_star[, b])

    # 2. Compute the "Meat" (Grouped scores)
    # Using rowsum is much faster than manual loops over clusters
    u_g <- rowsum(X * as.numeric(rb), gr)
    meat <- (t(u_g) %*% u_g) * df_c

    # 3. Assemble the Sandwich: Bread %*% Meat %*% Bread
    boot_vcov <- bread %*% meat %*% bread
    boot_se_matrix[, b] <- sqrt(diag(boot_vcov))
  }

  t_star <- (B_star - tb) / boot_se_matrix

  # for (i in 1:B) {
  #   w_raw <- sample(c(-1, 1), size = G, replace = TRUE)
  #   W <- w_raw[match(gr, clust_ids)]
  #   ystar <- yhat + re * W
  #   dat$ystar <- ystar
  #   if (!is.null(param)) {
  #     bhat <- obread %*% t(Xo) %*% ystar
  #     pred <- Xo %*% bhat
  #     resid2 <- as.numeric(ystar - pred)
  #     u_g <- rowsum(Xo * resid2, gr)
  #     meat <- (t(u_g) %*% u_g) * df_c
  #     boot_vcov <- obread %*% meat %*% obread
  #     boot_se <- sqrt(diag(boot_vcov))
  #     tsr[i] <- bhat[param, 1]/boot_se[param]
  #   }
  #   else {
  #     if (!cr3) {
  #       bhat <- bread %*% t(X) %*% ystar
  #       pred <- X %*% bhat
  #       resid2 <- as.numeric(ystar - pred)
  #       u_g <- rowsum(X * resid2, gr)
  #       meat <- (t(u_g) %*% u_g) * df_c
  #       boot_vcov <- bread %*% meat %*% bread
  #       boot_se <- sqrt(diag(boot_vcov))
  #       tstar[i, ] <- (bhat[, 1] - tb)/boot_se
  #     }
  #     else {
  #       m2 <- update(m1, ystar ~ .)
  #       bhat <- coef(m2)
  #       boot_se <- cr3(m2, dat = dat, cluster = cluster)
  #       tstar[i, ] <- (bhat - tb)/boot_se
  #     }
  #   }
  # }


  #########


  pv <- cil <- ciu <- numeric(k)

  pv <- numeric(k)
  for (j in 1:k) {
    # Using the (B+1) logic from earlier
    extreme_count <- sum(abs(t_star[j, ]) >= abs(tt[j]))
    pv[j] <- (extreme_count + 1) / (B + 1)
  }

  cil <- ciu <- numeric(k)
  for (j in 1:k) {
    # Get the critical t-values from the bootstrap distribution
    # Note: We use -pct[2] and -pct[1] because we are pivoting
    quants <- quantile(t_star[j, ], c(0.025, 0.975), na.rm = TRUE)
    cil[j] <- tb[j] - quants[2] * obsse[j]
    ciu[j] <- tb[j] - quants[1] * obsse[j]
  }

  # if (is.null(param)) {
  #   for (j in 1:k) {
  #     # pv[j] <- mean(abs(tstar[, j]) > abs(tt[j]))
  #     pct <- quantile(tstar[, j], c(0.025, 0.975), na.rm = TRUE)
  #     cil[j] <- tb[j] - pct[2] * obsse[j]
  #     ciu[j] <- tb[j] - pct[1] * obsse[j]
  #   }
  # }
  # else {
  #   pv <- mean(abs(tsr) > abs(tstat))
  #   pct <- quantile(tsr, c(0.025, 0.975), na.rm = TRUE)
  #   cil <- tb2 - pct[2] * obsse2[param]
  #   ciu <- tb2 - pct[1] * obsse2[param]
  # }
  if (is.null(param)) {
    res <- data.frame(term = names(tb), estimate = as.numeric(tb),
                      orig.std.error = as.numeric(obsse), wcb.p.value = pv,
                      conf.low = cil, conf.high = ciu)
  }
  else {
    res <- data.frame(term = param, estimate = tb2, wcb.p.value = pv,
                      conf.low = cil, conf.high = ciu)
  }
  rownames(res) <- NULL
  type = ifelse(cr3, "CR3", "CR1")
  type2 = ifelse(is.null(param), "WCU", "WCR")
  invisible(list(res = res, setype = type, wctype = type2))
}
