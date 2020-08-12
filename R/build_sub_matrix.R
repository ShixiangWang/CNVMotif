#' Build a Substitution Matrix
#'
#' @param simple_version if `TRUE`, just use segmental copy number value.
#' @param max_len_score the maximum score for segment length (should >=4).
#' The maximum score for copy number value is 6 (cannot be changed).
#'
#' @return a `list`.
#' @export
#'
#' @examples
#' sub_list <- build_sub_matrix()
#' sub_list2 <- build_sub_matrix(simple_version = TRUE)
#' @testexamples
#' expect_is(sub_list, "list")
#' expect_is(sub_list2, "list")
build_sub_matrix <- function(simple_version = FALSE, max_len_score = 4L) {
  stopifnot(max_len_score >= 4)

  if (isFALSE(simple_version)) {
    l <- 1:4
    v <- 0:5
    k <- LETTERS[1:24]
    map <- k
    names(map) <- vector_to_combination(1:4, 0:5)

    max_l <- max_len_score
    max_v <- length(v)
    pair_mat <- expand.grid(l, v, KEEP.OUT.ATTRS = FALSE) %>% as.matrix()

    score_mat <- pairScoreMatrix(pair_mat, pair_mat, max_l, max_v)
  } else {
    v <- 0:5
    k <- LETTERS[1:6]
    map <- k
    names(map) <- as.character(v)

    max_v <- length(v)
    pair_mat <- v %>% as.matrix()

    score_mat <- pairScoreSimpleMatrix(pair_mat, pair_mat, max_v)
  }

  rownames(score_mat) <- colnames(score_mat) <- k

  return(list(
    map = map,
    mat = score_mat
  ))
}
