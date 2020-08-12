#' Show Copy Number Sequence Logos
#' @inheritParams ggseqlogo2
#' @inheritParams build_sub_matrix
#' @param x a character vector of sequences or named list of sequences. All sequences must have same width.
#' @param recode if `TRUE`, it will use default indicator or specified indicator to show the letters in
#' the plot
#' @param indicator a named vector (like a dictory) to change letters one to one in the plot.
#' @return a `ggplot` object
#' @export
#' @examples
#' p1 <- show_seq_logo(sapply(split(LETTERS[1:24], 1:4), function(x) paste0(x, collapse = "")))
#' p1
#' p2 <- show_seq_logo(sapply(split(LETTERS[1:24], 1:4), function(x) paste0(x, collapse = "")),
#'   recode = TRUE
#' )
#' p2
#' p3 <- show_seq_logo(sapply(split(LETTERS[1:6], 1:2), function(x) paste0(x, collapse = "")),
#'   simple_version = TRUE
#' )
#' @testexamples
#' expect_is(p1, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(p3, "ggplot")
show_seq_logo <- function(x, method = c("prob", "bits"),
                          simple_version = FALSE,
                          ncol = NULL, nrow = NULL,
                          recode = FALSE, indicator = NULL, ...) {
  method <- match.arg(method)

  ## copy from sigminer utils.R
  reds <- sapply(list(c(252, 138, 106), c(241, 68, 50), c(188, 25, 26)),
    FUN = function(x) rgb2hex(x[1], x[2], x[3])
  ) %>% as.character()
  blues <- sapply(list(c(74, 152, 201), c(23, 100, 171)),
    FUN = function(x) rgb2hex(x[1], x[2], x[3])
  ) %>% as.character()

  if (simple_version) {
    cs <- ggseqlogo::make_col_scheme(
      chars = LETTERS[1:6],
      groups = c(
        rep("2 copy DEL", 1),
        rep("1 copy DEL", 1),
        rep("Normal", 1),
        rep("1 copy AMP", 1),
        rep("2 copy AMP", 1),
        rep("3+ copy AMP", 1)
      ),
      cols = c(
        rep("blue", 1),
        rep(blues[1], 1),
        rep("black", 1),
        rep(reds[1], 1),
        rep(reds[2], 1),
        rep(reds[3], 1)
      ),
      name = "Segment type"
    )

    ns <- LETTERS[1:6]
  } else {
    cs <- ggseqlogo::make_col_scheme(
      chars = LETTERS[1:24],
      groups = c(
        rep("2 copy DEL", 4),
        rep("1 copy DEL", 4),
        rep("Normal", 4),
        rep("1 copy AMP", 4),
        rep("2 copy AMP", 4),
        rep("3+ copy AMP", 4)
      ),
      cols = c(
        rep("blue", 4),
        rep(blues[1], 4),
        rep("black", 4),
        rep(reds[1], 4),
        rep(reds[2], 4),
        rep(reds[3], 4)
      ),
      name = "Segment type"
    )

    ns <- LETTERS[1:24]
  }

  if (recode) {
    if (is.null(indicator)) {
      if (isFALSE(simple_version)) {
        indicator <- rep(c("S", "M", "L", "E"), 6)
        names(indicator) <- LETTERS[1:24]
      }
    } else {
      if (is.null(names(indicator))) {
        stop("The indicator should have names to map.")
      }
    }
  } else {
    indicator <- NULL
  }

  ggseqlogo2(x,
    ncol = ncol,
    nrow = nrow,
    method = method,
    namespace = ns,
    col_scheme = cs,
    idor = indicator,
    ...
  )
}
