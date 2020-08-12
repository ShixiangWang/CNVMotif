#' Show Copy Number Sequence Shapes
#'
#' @inheritParams show_seq_logo
#' @inheritParams ggplot2::facet_wrap
#' @inheritParams build_sub_matrix
#' @param map default is `NULL`, a named string vector.
#' @param x_lab x lab.
#' @param y_lab y lab.
#' @param line_size_scale the scale size for line width.
#'
#' @return a `ggplot` object.
#' @export
#'
#' @examples
#' p <- show_seq_shape(c("ADGHK"))
#' p
#'
#' x <- list(a = c("ABCDE", "AXFDP"), b = c("KKDFH", "GKDFM"))
#' p2 <- show_seq_shape(x)
#' p2
#'
#' p3 <- show_seq_shape(c("ABCD"), simple_version = TRUE)
#' p3
#' @testexamples
#' expect_is(p, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(p3, "ggplot")
show_seq_shape <- function(x, map = NULL,
                           simple_version = FALSE,
                           line_size_scale = 3,
                           x_lab = ifelse(simple_version, "Assumed equal length", "Estimated segment length"),
                           y_lab = "Copy number",
                           nrow = NULL, ncol = NULL, scales = "free_x") {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required, please install it firstly!")
  }

  if (simple_version) {
    if (is.null(map)) {
      map <- 0:5
      names(map) <- LETTERS[1:6]
    }

    map_df <- data.frame(
      segVal = 0:5,
      stringsAsFactors = FALSE
    )
  } else {
    if (is.null(map)) {
      map <- vector_to_combination(1:4, 0:5)
      names(map) <- LETTERS[1:24]
    }

    map_df <- data.frame(
      lenVal = strsplit(map, split = "") %>% purrr::map_int(~ as.integer(.[1])),
      segVal = strsplit(map, split = "") %>% purrr::map_int(~ as.integer(.[2]))
    )
  }

  rownames(map_df) <- names(map)

  if (is.list(x)) {
    df <- purrr::map_df(x, seg_data, map_df = map_df, simple = simple_version, .id = "grp_id")
  } else {
    df <- seg_data(x, map_df, simple = simple_version)
  }

  p <- ggplot(df, aes_string(x = "x", y = "segVal", xend = "x_end", yend = "segVal")) +
    geom_segment(color = df$color, size = df$w * line_size_scale, alpha = df$w) +
    scale_y_continuous(breaks = 0:5, labels = c(0:4, "5+"), limits = c(0, 5)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = x_lab, y = y_lab) +
    ggseqlogo::theme_logo()

  if (is.list(x)) {
    p <- p + facet_wrap(~grp_id, nrow = nrow, ncol = ncol, scales = scales)
  }

  return(p)
}

seg_data <- function(x, map_df, simple = FALSE) {
  x_list <- strsplit(x, split = "")
  df <- purrr::map_df(x_list, function(x) {
    dplyr::tibble(
      pos = seq_along(x),
      seq = x
    )
  }, .id = "seq_id")
  # dplyr::count(.data$pos, .data$seq)
  df <- dplyr::bind_cols(df, map_df[df$seq, , drop = FALSE])

  df_freq <- df %>%
    dplyr::count(.data$pos, .data$seq) %>%
    dplyr::group_by(.data$pos) %>%
    dplyr::mutate(w = .data$n / sum(.data$n)) %>%
    dplyr::ungroup()

  if (simple) {
    df$lenVal <- 1L
  }
  df <- df %>%
    dplyr::group_by(.data$seq_id) %>%
    dplyr::mutate(
      x_end = cumsum(.data$lenVal),
      x = dplyr::lag(.data$x_end, default = 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(color = dplyr::case_when(
      .data$segVal > 2 ~ "red",
      .data$segVal < 2 ~ "blue",
      TRUE ~ "black"
    )) %>%
    dplyr::select(-.data$seq_id) %>%
    unique() %>%
    dplyr::left_join(df_freq, by = c("pos", "seq"))

  df
}
