#' robust_glmm (Robust standard errors for glmer/merMod objects)
#'
#' Create a matrix of replicate weights.
#'
#' @param model The glmer model object.
#'
#' @return Outputs a data.frame of results comparing model-based standard errors compared to robust standard errors. diag signifies how much larger the robust vs MBSE is.
#'
#' @examples
#' \dontrun{
#' library(lme4)
#' data(suspend, package = 'MLMusingR')
#' m1 <- glmer(sus ~ gpa + frpmp + male + (1|school), data = suspend, family = binomial)
#' robust_glmm(m1)
#' }
#'
#' @import merDeriv
#' @importFrom sandwich sandwich
#' @importFrom lme4 fixef
#' @importFrom stats pnorm symnum vcov
#' @export
robust_glmm <- function(model) {
  V <- sandwich::sandwich(model)
  p <- length(lme4::fixef(model))
  Vfx <- as.matrix(V)[1:p, 1:p, drop = FALSE]   # fixed-effect block
  tmp <- summary(model)$coef
  est <- tmp[,1]
  se <- sqrt(diag(Vfx))

  z <- est / se
  p_val <- 2 * stats::pnorm(-abs(z))
  mbse <- tmp[,2]

  stars <- stats::symnum(
    p_val,
    corr = FALSE,
    na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )

  data.frame(
    Estimate    = est,
    MBSE        = mbse,
    `Robust SE` = se,
    `z value`   = z,
    `Pr(>|z|)`  = p_val,
    ` `         = as.character(stars),
    diag        = se / mbse,
    check.names = FALSE
  )
}
