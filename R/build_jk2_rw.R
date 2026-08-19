#' build_jk2_rw (Create a matrix of replicate weights for JK2)
#'
#' Create a matrix of replicate weights.
#'
#' @param zone The replicate zone.
#' @param rep The replicate member (0 or 1)
#' @param weight The weight at level 1.
#'
#' @return Outputs a matrix that can be used in svrepdesign.
#'
#' @examples
#' \dontrun{
#' rw <- build_jk2_rw(df2$zone, df2$member, df2$totwgt)
#' }
#' @export
#'
build_jk2_rw <- function(zone, rep, weight) {
  zn      <- unique(zone)
  n_zones <- length(zn)
  n       <- length(weight)
  m       <- matrix(weight, nrow = n, ncol = n_zones)
  colnames(m) <- paste0("rw", seq_len(n_zones))
  for (h in seq_len(n_zones)) {
    in_zone <- zone == zn[h]
    m[in_zone & rep == 0, h] <- weight[in_zone & rep == 0] * 2
    m[in_zone & rep == 1, h] <- 0
  }
  m
}
