#' Coding Copy Number Segments with Letters.
#'
#' See [get_score_matrix()] for examples. See details for full description
#' of implementation.
#'
#' For complicated cases, letters are grouped as short (<50kb), mid (<500kb), long (<5Mb),
#' long (or extreme) long (>5Mb) segments.
#' - A B C D for copy number 0.
#' - E F G H for copy number 1.
#' - I J K L for copy number 2.
#' - M N O P for copy number 3.
#' - Q R S T for copy number 4.
#' - U V W X for copy number 5+.
#'
#' For simplified cases, letters are used to code only segment copy number value.
#' - A for copy number 0.
#' - B for copy number 1.
#' - C for copy number 2.
#' - D for copy number 3.
#' - E for copy number 4.
#' - F for copy number 5+.
#'
#' @param x a `CopyNumber` object or a `data.frame` with
#' at least 5 columns ("sample", "chromosome", "start", "end", "segVal").
#' @inheritParams build_sub_matrix
#'
#' @return a `list`.
#' @export
transform_seqs <- function(x, simple_version = FALSE, max_len_score = 4L) {
  if (inherits(x, "CopyNumber")) {
    x <- x@data
  } else {
    if (data.table::is.data.table(x)) {
      x <- data.table::copy(x)
    } else {
      x <- data.table::as.data.table(x)
    }
  }

  ## Make sure the order is right
  x <- x[order(sample, chromosome, start)]

  x[, segVal := ifelse(segVal > 5, 5, segVal) %>% as.integer()]
  if (isFALSE(simple_version)) {
    x[, lenVal := cut(end - start + 1L,
                      breaks = c(-Inf, 5e4, 5e5, 5e6, Inf),
                      labels = c("1", "2", "3", "4"),
                      right = FALSE
    ) %>% as.integer()]
  }

  x[, ID := paste(sample, chromosome, sep = ":")]

  sub_list <- build_sub_matrix(simple_version = simple_version, max_len_score = max_len_score)

  if (simple_version) {
    x$Seqs <- sub_list$map[as.character(x$segVal)]
  } else {
    x$Seqs <- sub_list$map[paste0(x$lenVal, x$segVal)]
  }

  return(
    list(
      dt = x,
      map = sub_list$map,
      mat = sub_list$mat
    )
  )
}
