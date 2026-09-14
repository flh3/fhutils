#' create_jk2_zones (Create a jacknife zones and replicates)
#' Create the zones and replicate member prior to making weights
#'
#' @param data dataset.
#' @param psu_id ID of the primary sampling unit.
#' @param sort_var The measure of size
#' @param stratum The stratum (if present)
#' @param seed For replicability.
#'
#' @return Outputs a new dataset with the zone and member number.
#' Dataset is sorted by MOS. PSUs are paired and then randomly split.
#' One PSU gets its weight multiplied by two, the other multiplied by zero.
#' If there is an odd number of units, the last PSU is split in two to form
#' a pseudo strata.
#'
#' @examples
#' \dontrun{
#' rw <- create_jk2_zones(df2, id, mos, stype, seed = 123)
#' }
#' @export
create_jk2_zones <- function(data, psu_id, sort_var, stratum = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Resolve column names (handles both bare symbols and strings)
  psu_col  <- if (is.character(substitute(psu_id))) psu_id else as.character(substitute(psu_id))
  sort_col <- if (is.character(substitute(sort_var))) sort_var else as.character(substitute(sort_var))

  # Handle optional stratum
  has_stratum <- !missing(stratum) && !is.null(substitute(stratum))
  if (has_stratum) {
    strat_col <- if (is.character(substitute(stratum))) stratum else as.character(substitute(stratum))
    strat_vec <- as.character(data[[strat_col]])
  } else {
    strat_col <- "..internal_stratum.."
    strat_vec <- rep("1", nrow(data))
    data[[strat_col]] <- strat_vec
  }

  # 1. Build initial unique PSU frame and sort
  psu_frame <- unique(data[, c(psu_col, strat_col, sort_col)])
  psu_frame <- psu_frame[order(psu_frame[[strat_col]], psu_frame[[sort_col]]), ]

  # Identify odd-count strata and extract their last PSU
  strat_split <- split(psu_frame, psu_frame[[strat_col]])
  odd_leftovers <- character(0)

  for (st in names(strat_split)) {
    sub_df <- strat_split[[st]]
    n_psus <- nrow(sub_df)
    if (n_psus %% 2 == 1) {
      odd_leftovers <- c(odd_leftovers, as.character(sub_df[[psu_col]][n_psus]))
    }
  }

  # 2. Split odd leftover PSUs into two quasi-PSUs
  data$quasi_psu <- as.character(data[[psu_col]])

  for (id in odd_leftovers) {
    row_idx <- which(data[[psu_col]] == id)
    n_rows  <- length(row_idx)

    if (n_rows > 1) {
      half_idx <- sample(row_idx, size = floor(n_rows / 2))
      data$quasi_psu[half_idx]                  <- paste0(id, "_a")
      data$quasi_psu[setdiff(row_idx, half_idx)] <- paste0(id, "_b")
    } else {
      data$quasi_psu[row_idx] <- paste0(id, "_a")
    }
  }

  # 3. Rebuild frame at quasi-PSU level and assign JK2 zones + rep flags
  quasi_frame <- unique(data[, c("quasi_psu", psu_col, strat_col, sort_col)])
  quasi_frame <- quasi_frame[order(quasi_frame[[strat_col]],
                                   quasi_frame[[sort_col]],
                                   quasi_frame$quasi_psu), ]

  # Pair units within stratum
  quasi_split <- split(quasi_frame, quasi_frame[[strat_col]])

  paired_list <- lapply(quasi_split, function(sub_df) {
    n_units <- nrow(sub_df)
    pair_id <- ceiling(seq_len(n_units) / 2)

    # Format zone prefix nicely if unstratified
    if (has_stratum) {
      sub_df$zone <- paste(sub_df[[strat_col]], pair_id, sep = "_")
    } else {
      sub_df$zone <- paste0("zone_", pair_id)
    }

    # Assign 0 or 1 indicator within each zone
    zone_counts <- table(sub_df$zone)
    rep_vals <- unlist(lapply(zone_counts, function(k) seq_len(k) - 1))
    sub_df$rep <- as.vector(rep_vals)

    sub_df[, c("quasi_psu", "zone", "rep")]
  })

  lookup_table <- do.call(rbind, paired_list)
  rownames(lookup_table) <- NULL

  # 4. Clean up temporary column if unstratified
  if (!has_stratum) {
    data[[strat_col]] <- NULL
  }

  # Merge zone and rep indicator back onto the original data
  merge(data, lookup_table, by = "quasi_psu", all.x = TRUE, sort = FALSE)
}
