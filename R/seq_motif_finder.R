#' Coding Copy Number Segments with Letters.
#'
#' See [sh_get_score_matrix()] for examples.
#' Letters are grouped as short (<50kb), mid (<500kb), long (<5Mb),
#' long (or extreme) long (>5Mb) segments.
#' - A B C D for copy number 0.
#' - E F G H for copy number 1.
#' - I J K L for copy number 2.
#' - M N O P for copy number 3.
#' - Q R S T for copy number 4.
#' - U V W X for copy number 5+.
#'
#' @param x a `CopyNumber` object or a `data.frame` with
#' at least 5 columns ("sample", "chromosome", "start", "end", "segVal").
#' @inheritParams sh_build_sub_matrix
#'
#' @return a `list`.
#' @export
sh_coding_segs <- function(x, max_len_score = 8L) {
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

  x[, `:=`(
    lenVal = end - start + 1L,
    segVal = ifelse(segVal > 5, 5, segVal) %>% as.integer() ## Set max value
  )]
  x[, lenVal := cut(lenVal,
    breaks = c(-Inf, 5e4, 5e5, 5e6, Inf),
    labels = c("1", "2", "3", "4"),
    right = FALSE
  ) %>% as.integer()]

  x[, ID := paste(sample, chromosome, sep = ":")]

  sub_list <- sh_build_sub_matrix(max_len_score = max_len_score)
  x$Seqs <- sub_list$map[paste0(x$lenVal, x$segVal)]
  x$Seqs

  return(
    list(
      dt = x,
      map = sub_list$map,
      mat = sub_list$mat
    )
  )
}

#' Build a Substitution Matrix
#'
#' @param max_len_score the maximum score for segment length (should >=4).
#' Default is 8 for balancing the weights between segment length and copy number
#' value. The maximum score for copy number value is 6.
#'
#' @return a `list`.
#' @export
#'
#' @examples
#' sub_list <- sh_build_sub_matrix()
#' @testexamples
#' expect_is(sub_list, "list")
sh_build_sub_matrix <- function(max_len_score = 8L) {
  stopifnot(max_len_score >= 4)

  l <- 1:4
  v <- 0:5
  k <- LETTERS[1:24]
  map <- k
  names(map) <- vector_to_combination(1:4, 0:5)

  max_l <- max_len_score
  max_v <- length(v)
  pair_mat <- expand.grid(l, v, KEEP.OUT.ATTRS = FALSE) %>% as.matrix()

  score_mat <- pairScoreMatrix(pair_mat, pair_mat, max_l, max_v)
  rownames(score_mat) <- colnames(score_mat) <- k

  return(list(
    map = map,
    mat = score_mat
  ))
}

collapse_shift_seqs <- function(x, len = 5L, step = 1L) {
  if (length(x) <= len) {
    return(paste(x, collapse = ""))
  } else {
    y <- c()
    i <- 0L
    while (i <= (length(x) - len)) {
      y <- c(y, paste(x[(1 + i):(len + i)], collapse = ""))
      i <- i + step
    }
    return(y)
  }
}

#' Extract Pasted Sequences from Each Chromosome
#'
#' See [sh_get_score_matrix()] for examples.
#' The result sequences are unique and sorted.
#'
#' @param dt a `data.table` from [sh_coding_segs].
#' @param len cut length.
#' @param step step size to move on each chromosome sequence.
#' @param return_dt if `TRUE`, just return a `data.table` containing
#' mutated `Seqs` column.
#'
#' @return a `list`.
#' @export
sh_extract_seqs <- function(dt, len = 5L, step = 2L, return_dt = FALSE) {
  stopifnot(data.table::is.data.table(dt))

  dt <- dt[, c("ID", "Seqs")]
  dt <- dt[, list(Seqs = collapse_shift_seqs(Seqs, len = len, step = step)),
    by = "ID"
  ]

  if (return_dt) {
    return(dt)
  }

  all_seqs <- unique(dt$Seqs)
  keep <- nchar(all_seqs) >= len

  list(
    keep = sort(all_seqs[keep]),
    drop = sort(all_seqs[!keep])
  )
}

#' Get Copy Number Sequence Similarity or Distance Matrix
#'
#' @param x a coding copy number sequence (valid letters are A to X).
#' @param sub_mat a substitution matrix, each element indicates a score to plus.
#' See [sh_build_sub_matrix()].
#' @param block_size a block size to aggregrate, this is designed for big data, it means
#' results from adjacent sequences will be aggregrate by means to reduce the size of result
#' matrix.
#' @param verbose if `TRUE`, print extra message, note it will slower the computation.
#' @param cores computer cores, default is `1`, note it is super fast already, set more
#' cores typically do not speed up the computation.
#'
#' @return a score matrix.
#' @export
#'
#' @examples
#' load(system.file("extdata", "toy_segTab.RData",
#'   package = "sigminer.helper", mustWork = TRUE
#' ))
#' x <- sh_coding_segs(segTabs)
#' x
#' seqs <- sh_extract_seqs(x$dt)
#' seqs
#' mat <- sh_get_score_matrix(seqs$keep, x$mat, verbose = TRUE)
#' mat[1:5, 1:5]
#'
#' mat2 <- sh_get_score_matrix2(seqs$keep, x$mat)
#' identical(mat, mat2)
#'
#' mat_b <- sh_get_score_matrix(seqs$keep, x$mat, block_size = 2L)
#' ## block1 represents the first 2 sequences
#' ## block2 represents the 3rd, 4th sequences
#' ## ...
#' mat_b[1:5, 1:5]
#'
#' if (require("doParallel")) {
#'   mock_seqs <- sapply(1:10000, function(x) {
#'     paste(sample(LETTERS[1:24], 5, replace = TRUE), collapse = "")
#'   })
#'
#'   system.time(
#'     y1 <- sh_get_score_matrix(mock_seqs, x$mat, cores = 1)
#'   )
#'
#'   system.time(
#'     y2 <- sh_get_score_matrix(mock_seqs, x$mat, cores = 2)
#'   )
#'
#'   all.equal(y1, y2)
#' }
#' @testexamples
#' expect_is(x, "list")
#' expect_is(seqs, "list")
#' expect_is(mat, "matrix")
#' expect_equal(mat, mat2)
#' if (require("doParallel")) {
#'   expect_equal(y1, y2)
#' }
sh_get_score_matrix <- function(x, sub_mat, block_size = NULL, verbose = FALSE, cores = 1L) {
  stopifnot(is.numeric(cores))

  if (anyNA(sub_mat)) {
    stop("Input substitution matrix cannot contain 'NA' values!")
  }

  map <- seq_len(24L)
  names(map) <- LETTERS[map]
  map <- map - 1L # to 0 based index

  ## Checking input
  if (any(grepl("[^A-X]", x, ignore.case = FALSE))) {
    stop("The input sequences should contain only A->X, any other letters are invalid.")
  }

  m <- matrix(NA_integer_, ncol = length(x), nrow = nchar(x[1]))

  for (i in seq_len(ncol(m))) {
    s <- unlist(strsplit(x[i], split = ""))
    m[, i] <- map[s] %>% as.integer()
  }
  m <- t(m)

  if (!is.null(block_size)) {
    stopifnot(block_size > 1)
  } else {
    block_size <- 1
  }

  if (cores == 1) {
    y <- getScoreMatrix(m, sub_mat, block_size, verbose)

    if (block_size == 1) {
      colnames(y) <- rownames(y) <- x
    } else {
      colnames(y) <- rownames(y) <- paste0("block", seq_len(nrow(y)))
    }
  } else {
    if (block_size > 1) {
      stop("In parallel mode, 'block_size' can only be one!")
    }

    if (nrow(m) < 10000) {
      warning("For data <10000, set cores > 1 is not recommended.", immediate. = TRUE)
    }

    if (cores <= 0 | cores > parallel::detectCores()) {
      cores <- parallel::detectCores() %>% as.integer()
    }

    ngrp <- ceiling(nrow(m) / 1000)
    grp_list <- chunk2(seq_len(nrow(m)), ngrp)

    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'doParallel' is required to go through this parallel method.")
    }

    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required to go through this parallel method.")
    }

    if (Sys.info()[["sysname"]] == "Windows") {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
    } else {
      doParallel::registerDoParallel(cores = cores)
    }

    y <- foreach(
      i = seq_along(grp_list),
      .combine = "cbind",
      .packages = "sigminer.helper",
      .export = c("m", "grp_list", "sub_mat", "verbose", "getScoreMatrixRect"),
      .verbose = FALSE
    ) %dopar% {
      getScoreMatrixRect(m, m[grp_list[[i]], ], sub_mat, verbose)
    }

    # for (i in seq_along(grp_list)) {
    #   y[, grp_list[[i]]] <- getScoreMatrixRect(m, m[grp_list[[i]], ], sub_mat, verbose)
    # }

    colnames(y) <- rownames(y) <- x
  }

  return(y)
}

score_pairwise_strings <- function(x, y, sub_mat) {
  ## For length-1 string "ABCD"
  ## or vector c("A", "B", "C", "D")
  if (length(x) == 1) {
    x <- strsplit(x, "")[[1]]
  }
  if (length(y) == 1) {
    y <- strsplit(y, "")[[1]]
  }
  sub_mat[x, y] %>%
    diag() %>%
    sum()
}

#' `sh_get_score_matrix2()` is a variant version of `sh_get_score_matrix()`. Normally, it worker worse than `sh_get_score_matrix()`.
#' @rdname sh_get_score_matrix
#' @param method a method for getting (storing) the results.
#' @export
sh_get_score_matrix2 <- function(x, sub_mat, method = c("base", "ff", "bigmemory"), verbose = FALSE) {
  method <- match.arg(method)
  n <- length(x)

  if (method == "base") {
    mat <- matrix(NA_integer_, nrow = n, ncol = n)
  } else if (method == "ff") {
    mat <- ff::ff(NA_integer_,
      dim = c(n, n), vmode = "byte"
    ) ## Byte from -128 ~ 127
  } else {
    options(bigmemory.allow.dimnames = TRUE, bigmemory.typecast.warning = FALSE)
    mat <- bigmemory::big.matrix(n, n, type = "integer")
  }
  # Matrix column is faster than row

  i <- j <- 1
  for (i in seq_len(n)) {
    if (verbose) message("Handling sequence: ", x[i])
    j_vals <- vector(mode = "integer", length = i)
    for (j in seq_len(i)) {
      j_vals[j] <- score_pairwise_strings(x[i], x[j], sub_mat = sub_mat) %>% as.integer()
    }
    mat[seq_along(j_vals), i] <- j_vals
  }

  mat <- mat[]
  ## NOTE the t() operation
  ## cannot just assign upper to lower triangle matrix
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  rownames(mat) <- colnames(mat) <- x
  return(mat)
}

show_seq_shape <- function(x, map = NULL, x_lab = "Estimated segment length", y_lab = "Copy number") {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required, please install it firstly!")
  }

  if (is.null(map)) {
    map <- LETTERS[1:24]
    names(map) <- vector_to_combination(1:4, 0:5)
  }

  map_df <- data.frame(
    lenVal = strsplit(names(map), split = "") %>% purrr::map_int(~ as.integer(.[1])),
    segVal = strsplit(names(map), split = "") %>% purrr::map_int(~ as.integer(.[2])),
    stringsAsFactors = FALSE
  )
  rownames(map_df) <- as.character(map)

  ## test data
  df <- map_df[unlist(strsplit(x, split = "")), ]
  df$x_end <- cumsum(df$lenVal)
  df$x <- dplyr::lag(df$x_end, default = 0)
  df$color <- ifelse(df$segVal > 2, "red",
    ifelse(df$segVal < 2, "blue",
      "black"
    )
  )

  ggplot(df, aes_string(x = "x", y = "segVal", xend = "x_end", yend = "segVal")) +
    geom_segment(color = df$color) +
    scale_y_continuous(breaks = 0:5, labels = c(0:4, "5+"), limits = c(0, 5)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = x_lab, y = y_lab) +
    cowplot::theme_cowplot()
}

#' Show Copy Number Sequence Logos
#' @inheritParams ggseqlogo2
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
#' @testexamples
#' expect_is(p1, "ggplot")
#' expect_is(p2, "ggplot")
show_seq_logo <- function(x, method = c("prob", "bits"), ncol = NULL, nrow = NULL,
                          recode = FALSE, indicator = NULL, ...) {
  method <- match.arg(method)

  ## copy from sigminer utils.R
  reds <- sapply(list(c(252, 138, 106), c(241, 68, 50), c(188, 25, 26)),
    FUN = function(x) rgb2hex(x[1], x[2], x[3])
  ) %>% as.character()
  blues <- sapply(list(c(74, 152, 201), c(23, 100, 171)),
    FUN = function(x) rgb2hex(x[1], x[2], x[3])
  ) %>% as.character()

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

  if (recode) {
    if (is.null(indicator)) {
      indicator <- rep(as.character(1:4), 6)
      names(indicator) <- LETTERS[1:24]
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
    namespace = LETTERS[1:24],
    col_scheme = cs,
    idor = indicator,
    ...
  )
}
