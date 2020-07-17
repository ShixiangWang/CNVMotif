vector_to_combination <- function(..., c_string = "") {
  expand.grid(
    ...,
    stringsAsFactors = FALSE
  ) %>%
    apply(1, paste0, collapse = c_string) %>%
    unique()
}

# From https://gist.github.com/mbannert/e9fcfa86de3b06068c83
rgb2hex <- function(r, g, b) grDevices::rgb(r, g, b, maxColorValue = 255)
col2hex <- function(col, alpha) grDevices::rgb(t(grDevices::col2rgb(col)), alpha = alpha, maxColorValue = 255)

# https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))
