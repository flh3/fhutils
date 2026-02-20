#' ctab (Create cross tables)
#'
#' Creates cross tables showing the frequencies and percentages.
#'
#' @param data The dataframe.
#' @param row_var Variable on the y axis.
#' @param col_var Variable on the x axis.
#' @param percent_type Show percent by "row", "col", or "all".
#' @param show_na Whether to show missing values or not.
#'
#' @return Outputs a 'tabyl' / 'data.frame' object.
#'
#' @examples
#' ctab(mtcars, am, vs)
#'
#' @importFrom rlang .data
#' @importFrom janitor tabyl adorn_totals adorn_percentages adorn_ns adorn_title adorn_pct_formatting
#' @export
#'
ctab <- function(data, row_var, col_var, percent_type = "row",
                 show_na = TRUE) {

  # require(janitor)
  # require(dplyr)

  data |>
    tabyl({{ row_var }}, {{ col_var }}, show_na = show_na) |>
    adorn_totals(c("row", "col")) |>
    adorn_percentages(percent_type) |>
    adorn_pct_formatting(digits = 1) |>
    adorn_ns(position = "front") |>
    adorn_title("combined")

}
