#' Split Cluster Sequence into List
#'
#' @param x a named integer vector from `hclust` etc.
#' @param s default is `NULL`, when the `x` is block cluster sequence, set this
#' to a sequence vector.
#' @param block_size block size used to split, only used when `s` is not `NULL`.
#'
#' @return a `list`.
#' @export
cluster_split <- function(x, s = NULL, block_size = 10) {
  if (as.character(startsWith(x[1]), "block")) {
    cluster_blocks <- x
    dt <- dplyr::tibble(
      block = names(cluster_blocks),
      cluster = as.character(cluster_blocks),
      seq = chunk2(s, ceiling(length(s) / block_size))
    ) %>%
      tidyr::unnest("seq") %>%
      data.table::as.data.table()

    split(dt$seq, dt$cluster)
  } else {
    split(names(x), as.integer(x))
  }
}
