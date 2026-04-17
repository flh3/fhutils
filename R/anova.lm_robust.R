#' anova (model comparison)
#'
#' @description Computes a deviance table for model comparisons using lm_robust objects.
#'
#' @param ... lm_robust objects containing the result from the lm_robust function.
#'
#' @return Outputs the model comparison table (as a data.frame). If not statistically significant,
#' the simpler model may be preferred.
#'
#' @section Warning:
#' This function requires that the model objects be entered from the simplest
#' to the more complex model (reduced model first, then full model). The function
#' requires that the models be nested within each other.
#'
#' @examples
#' \dontrun{
#' library(estimatr)
#' m1 <- lm_robust(mpg ~ wt, data = mtcars)
#' m2 <- lm_robust(mpg ~ wt + am, data = mtcars)
#' m3 <- lm_robust(mpg ~ wt + am + vs, data = mtcars)
#'
#' anova(m1, m2, m3)
#' }
#'
#' @importFrom stats pf symnum
#' @export
#' @method anova lm_robust
anova.lm_robust <- function(...) {

  mods <- list(...)
  n_mods <- length(mods)

  if (n_mods < 2) stop("Provide at least two models for comparison.")

  # Check sample sizes across all models
  n_obs <- sapply(mods, nobs)
  if (length(unique(n_obs)) > 1) stop("Models must have the same sample size.")

  # Initialize results table
  res <- data.frame(
    df_reg = sapply(mods, function(x) length(x$term)),
    df_res = sapply(mods, function(x) nobs(x) - length(x$term)),
    R2     = sapply(mods, function(x) x$r.squared)
  )

  res$df_diff <- c(NA, diff(res$df_reg))
  res$R2_diff <- c(NA, diff(res$R2))

  if (any(res$df_diff < 0, na.rm = T)) stop('Enter models from simpler (reduced) to more complex (full).')

  # Calculate F-stats and p-values
  # F = ((R2_u - R2_r) / df_diff) / ((1 - R2_u) / df_res_u)
  res$Fstat <- round((res$R2_diff / res$df_diff) / ((1 - res$R2) / res$df_res), 3)
  res$p     <- round(pf(res$Fstat, df1 = res$df_diff, df2 = res$df_res, lower.tail = FALSE), 4)

  res$" " = as.character(symnum(res$p,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
         symbols   = c("***", "**", "*", ".", " "),
         abbr.colnames = FALSE,
         na = FALSE))

  class(res) <- c("anova_robust2", "data.frame")
  return(res)
}


