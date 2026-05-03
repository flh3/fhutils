#' pool_ps_subclass (Pooling function for propensity scores using sub classifications)
#'
#' Computes the estimand of interest using subclassifications with propensity scores.
#' A much simpler way would just be to use the weights produced by the MatchIt package when using the subclass option.
#'
#' @param data The dataset to use.
#' @param outcome The outcome of interest.
#' @param treatment The treatment variable.
#' @param subclass The subclassification variable.
#'
#' @return Outputs the treatment effect.
#' @importFrom stats var
#' @export
pool_ps_subclass <- function(data, outcome, treatment, subclass) {
  # Split the data into a list of data frames, one for each stratum
  strata_data <- split(data, data[[subclass]])

  # Initialize vectors to store the stratum-specific values
  k <- length(strata_data)
  tau_k <- numeric(k)
  var_tau_k <- numeric(k)
  N_k <- numeric(k)
  N1_k <- numeric(k)


  # Loop through each stratum to calculate local effects and variances
  for (i in seq_along(strata_data)) {
    df <- strata_data[[i]]
    if (is.factor(df[[treatment]])) df[[treatment]] <- as.numeric(df[[treatment]]) - 1
    # Extract outcomes for treated (1) and control (0) units
    y1 <- df[[outcome]][df[[treatment]] == 1]
    y0 <- df[[outcome]][df[[treatment]] == 0]

    # Calculate sample sizes
    n1 <- length(y1)
    n0 <- length(y0)
    N_k[i] <- n1 + n0
    N1_k[i] <- n1

    # Calculate means
    m1 <- mean(y1, na.rm = TRUE)
    m0 <- mean(y0, na.rm = TRUE)

    # Calculate variances (requires at least 2 observations per group)
    v1 <- ifelse(n1 > 1, var(y1, na.rm = TRUE), 0)
    v0 <- ifelse(n0 > 1, var(y0, na.rm = TRUE), 0)

    # Stratum-specific treatment effect (tau_k) and its variance
    tau_k[i] <- m1 - m0
    var_tau_k[i] <- (v1 / n1) + (v0 / n0)
  }

  # Calculate totals for weighting
  N_total <- sum(N_k)
  N1_total <- sum(N1_k)

  # ---------------------------------------------------------
  # 1. Pool for Average Treatment Effect (ATE)
  # ---------------------------------------------------------
  ate_weights <- N_k / N_total
  ate <- sum(ate_weights * tau_k)
  ate_var <- sum((ate_weights^2) * var_tau_k)
  ate_se <- sqrt(ate_var)

  # ---------------------------------------------------------
  # 2. Pool for Average Treatment Effect on the Treated (ATT)
  # ---------------------------------------------------------
  att_weights <- N1_k / N1_total
  att <- sum(att_weights * tau_k)
  att_var <- sum((att_weights^2) * var_tau_k)
  att_se <- sqrt(att_var)

  # ---------------------------------------------------------
  # Assemble the final results table
  # ---------------------------------------------------------
  results <- data.frame(
    Estimand = c("ATE", "ATT"),
    Estimate = c(ate, att),
    Std_Error = c(ate_se, att_se),
    CI_Lower_95 = c(ate - 1.96 * ate_se, att - 1.96 * att_se),
    CI_Upper_95 = c(ate + 1.96 * ate_se, att + 1.96 * att_se)
  )

  return(results)
}
