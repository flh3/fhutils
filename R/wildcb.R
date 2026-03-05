#' wildcb
#'
#' Fits a wild cluster bootstrapped for an lm model.
#'
#' @param fml A regression formula.
#' @param dat Original dataset.
#' @param cluster The clustering variable (in quotes).
#' @param B Number of bootstrap replicates.
#' @param seed A seed for replicability.
#' @param cr3 Whether to use the CR3 standard error correction. Only available for unrestricted models.
#' @param param Specify the coefficient of interest for a restricted model.
#'
#' @return Regression results.
#'
#' @examples
#' g1 <- wildcb(mpg ~ wt + am, cluster = 'cyl', dat = mtcars)
#
#'
#' @importFrom stats coef lm model.frame quantile update as.formula
#' @export
#'


wildcb <- function(fml, cluster, dat, B = 999, seed = 0, cr3 = FALSE,
                   param = NULL) {

  if (!is.null(param)){ #WCR
    m1  <- lm(fml, data = dat)

    n <- nobs(m1)
    Xo <- model.matrix(m1)
    # Xt <- t(X)
    obread <- vcov(m1) / (summary(m1)$sigma^2)

    if(!cr3) { #cr1
      obsse2 <- cr1(m1, cluster = cluster, dat = dat)
    } else { #cr3
      obsse2 <- cr3(m1, cluster = cluster, dat = dat)
    }

    tb2 = coef(m1)[param]
    tstat = tb2 / obsse2[param]

    drop_var <- param
    fml <- update(fml, as.formula(paste("~ . -", drop_var)))
  }

  # 1. Initial Model and Parameters
  m1  <- lm(fml, data = dat)
  tb  <- coef(m1)
  gr  <- dat[[cluster]]
  n <- nobs(m1) #how many
  y <- model.frame(m1)[[1]] #outcome

  if(!cr3) { #cr1
    obsse <- cr1(m1, cluster = cluster, dat = dat)
  } else { #cr3
    obsse <- cr3(m1, cluster = cluster, dat = dat)
  }

  tt  <- tb / obsse #t value

  yhat <- m1$fitted.values #fitted(m1)
  re  <- m1$residuals #residuals(m1)
  clust_ids <- unique(gr)
  G <- length(clust_ids)
  k  <- length(tb)

  tstar <- matrix(NA, ncol = k, nrow = B) #for WCU
  tsr <- numeric() #for WCR
  X <- model.matrix(m1)
  bread <- vcov(m1) / (summary(m1)$sigma^2) # = to (X'X)^-1
  if (!is.null(param)) k = k + 1
  df_c <- ((n - 1) / (n - k)) * (G / (G - 1)) #adjustment

  # 2. Bootstrap Loop
  if(seed != 0) set.seed(seed)
  for (i in 1:B) {

    w_raw <- sample(c(-1, 1), size = G, replace = TRUE)
    W <- w_raw[match(gr, clust_ids)]

    ystar <- yhat + re * W
    dat$ystar <- ystar #for lm.fit to work
    # m2 <- update(m1, ystar ~ .)

    if (!is.null(param)){ #WCR, only CR1
      # m2 <- lm.fit(Xo, ystar)

      bhat <- obread %*% t(Xo) %*% ystar
      pred <- Xo %*% bhat
      resid2 <- as.numeric(ystar - pred)
      # u_g <- rowsum(Xo * m2$residuals, gr)
      u_g <- rowsum(Xo * resid2, gr)
      meat <- (t(u_g) %*% u_g) * df_c
      boot_vcov <- obread %*% meat %*% obread
      boot_se <- sqrt(diag(boot_vcov))

      # tsr[i] <- m2$coefficients[param] / boot_se[param]
      tsr[i] <- bhat[param, 1] / boot_se[param]

    } else { #WCU

      if (!cr3){ #cr1

        bhat <- bread %*% t(X) %*% ystar
        pred <- X %*% bhat
        resid2 <- as.numeric(ystar - pred)
        # m2 <- lm.fit(X, ystar) #faster than lm
        # u_g <- rowsum(X * m2$residuals, gr)
        u_g <- rowsum(X * resid2, gr)
        meat <- (t(u_g) %*% u_g) * df_c
        boot_vcov <- bread %*% meat %*% bread
        boot_se <- sqrt(diag(boot_vcov))
      } else { #cr3
        m2 <- update(m1, ystar ~ .)
        boot_se <- cr3(m2, dat = dat, cluster = cluster)
      }
      tstar[i, ] <- (bhat[, 1] - tb) / boot_se
    }

  }


  # 3. Calculations
  pv <- cil <- ciu <- numeric(k)

  if (is.null(param)){ # unrestricted
    for (j in 1:k) {
      pv[j]  <- mean(abs(tstar[, j]) > abs(tt[j]))
      pct    <- quantile(tstar[, j], c(0.025, 0.975), na.rm = TRUE)
      cil[j] <- tb[j] - pct[2] * obsse[j] #eq. 8 MacKinnon
      ciu[j] <- tb[j] - pct[1] * obsse[j]
    }
  } else { # restricted

    pv  <- mean(abs(tsr) > abs(tstat))
    pct <- quantile(tsr, c(0.025, 0.975), na.rm = TRUE)
    cil <- tb2 - pct[2] * obsse2[param] #eq. 8 MacKinnon
    ciu <- tb2 - pct[1] * obsse2[param]

  }


  # 4. Results Table
  if (is.null(param)){
    res <- data.frame(
      term = names(tb),
      estimate = as.numeric(tb),
      orig.std.error = as.numeric(obsse),
      wcb.p.value = pv,
      conf.low = cil,
      conf.high = ciu
    )
  } else {
    res <- data.frame(
      term = param,
      estimate = tb2,
      wcb.p.value = pv,
      conf.low = cil,
      conf.high = ciu
    )
  }

  rownames(res) <- NULL
  type = ifelse(cr3, "CR3", 'CR1')
  type = ifelse(is.null(param), type, "CR1") #only available
  type2 = ifelse(is.null(param), 'WCU', 'WCR')
  # cat("Standard Error:", type, '\n')
  # cat("WCB Type:", type2, '\n')
  # print(res)
  invisible(list(res = res, setype = type, wctype = type2))

}
