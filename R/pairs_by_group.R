#' Generate All Pairs Within Each Group
#'
#' For each group in a data frame, generate all \eqn{C(n, 2)} pairwise
#' combinations of values in a specified ID column. Useful for KaKs
#' calculation between gene pairs within clades, co-evolution analysis,
#' or within-group network construction.
#'
#' @param df A data frame.
#' @param group_col Column to group by (bare name or character string).
#' @param id_col Column containing IDs to pair (bare name or character string).
#' @param out_names Character vector of length 2 naming the two output pair
#'   columns (default: \code{c("Gene1", "Gene2")}).
#'
#' @return A tibble with columns: \code{group_col}, \code{out_names[1]},
#'   \code{out_names[2]}. Groups with fewer than 2 IDs are silently dropped.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   type = c("A", "A", "A", "B", "B"),
#'   gene = c("g1", "g2", "g3", "g4", "g5")
#' )
#' pairs_by_group(df, type, gene)
#' # # A tibble: 4 x 3
#' #   type  Gene1 Gene2
#' #   <chr> <chr> <chr>
#' # 1 A     g1    g2
#' # 2 A     g1    g3
#' # 3 A     g2    g3
#' # 4 B     g4    g5
#'
#' # Original inline pattern this replaces:
#' # df %>% group_by(type) %>% group_modify(~ { ... combn(ids, 2) ... })
#' }
#'
pairs_by_group <- function(df, group_col, id_col,
                           out_names = c("Gene1", "Gene2")) {

  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  if (!is.character(out_names) || length(out_names) != 2 ||
      any(is.na(out_names)) || any(!nzchar(out_names))) {
    stop("'out_names' must be a character vector of length 2 with non-empty strings")
  }

  id_name <- rlang::as_name(rlang::enquo(id_col))

  df %>%
    dplyr::group_by({{ group_col }}) %>%
    dplyr::group_modify(~ {
      ids <- .x[[id_name]]
      if (length(ids) < 2) {
        return(stats::setNames(
          data.frame(character(0), character(0), stringsAsFactors = FALSE),
          out_names
        ))
      }
      if (length(ids) == 2) {
        return(stats::setNames(
          data.frame(t(ids), stringsAsFactors = FALSE),
          out_names
        ))
      }
      combos <- utils::combn(ids, 2, simplify = FALSE)
      stats::setNames(
        as.data.frame(do.call(rbind, combos), stringsAsFactors = FALSE),
        out_names
      )
    }) %>%
    dplyr::ungroup()
}
