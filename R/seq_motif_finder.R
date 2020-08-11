#' Coding Copy Number Segments with Letters.
#'
#' See [sh_get_score_matrix()] for examples. See details for full description
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
#' @inheritParams sh_build_sub_matrix
#'
#' @return a `list`.
#' @export
sh_coding_segs <- function(x, simple_version = FALSE, max_len_score = 4L) {
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

  sub_list <- sh_build_sub_matrix(simple_version = simple_version, max_len_score = max_len_score)

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
#' sub_list <- sh_build_sub_matrix()
#' sub_list2 <- sh_build_sub_matrix(simple_version = TRUE)
#' @testexamples
#' expect_is(sub_list, "list")
#' expect_is(sub_list2, "list")
sh_build_sub_matrix <- function(simple_version = FALSE, max_len_score = 4L) {
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

collapse_shift_seqs2 <- function(x, too_large, len = 5L, step = 1L) {
  z <- seq_along(x)[!too_large]
  if (length(z) == 0) {
    return(NULL)
  } else if (length(z) == 1) {
    return(x[z])
  } else {
    z_list <- split(x[z], findInterval(z, z[diff(z) > 1] + 2L))
    return(sapply(z_list, collapse_shift_seqs, len = len, step = step) %>% unlist() %>% as.character())
  }
}

collapse_local_seqs <- function(x, too_large, segsize, cutoff = 1e7) {
  z <- seq_along(x)[!too_large]
  zL <- segsize[!too_large]
  if (length(z) == 0) {
    return(NULL)
  } else if (length(z) == 1) {
    return(x[z])
  } else {
    iv <- findInterval(z, z[diff(z) > 1] + 2L)
    z_list <- split(x[z], iv)
    zL_list <- split(zL, iv)
    # Loop both segment copy number value and its length
    purrr::map2(z_list, zL_list, function(x, y, cutoff) {
      z <- getLocalSubstr(x, y, cutoff)
      z <- setdiff(z, "")  # Remove blank string
      z
    }, cutoff = cutoff) %>% purrr::flatten_chr()
    #return(sapply(z_list, paste, collapse = "") %>% unlist() %>% as.character())
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
#' @param local_cutoff any segment with length greater than this cutoff will be filtered out and
#' used as cutpoint, default is `10Mb`.
#' @param flexible_approach if `TRUE`, extract flexible-size sequences between segments
#'  with size less than specified cutoff. So the arguments `len` and `step` are ignored.
#' @param return_dt if `TRUE`, just return a `data.table` containing
#' mutated `Seqs` column.
#'
#' @return a `list`.
#' @export
sh_extract_seqs <- function(dt, len = 5L, step = 1L, local_cutoff = 1e7,
                            flexible_approach = FALSE,
                            return_dt = FALSE) {
  stopifnot(data.table::is.data.table(dt))

  if (flexible_approach) {
    message("Task: extract flexible-size sequences between segments with total size less than specified cutoff.")
    dt$segsize <- dt$end - dt$start + 1L
    dt$too_large <- dt$segsize >= local_cutoff
    message("Total segments is ", nrow(dt), " and ", sum(dt$too_large), " of them with length >=", local_cutoff)
    message("Fraction: ", round(sum(dt$too_large) / nrow(dt), digits = 3))
    dt <- dt[, c("ID", "Seqs", "too_large", "segsize")]

    dt <- dt[, list(Seqs = collapse_local_seqs(Seqs, too_large, segsize, cutoff = local_cutoff)), by = "ID"]

    if (return_dt) {
      return(dt)
    } else {
      return(sort(unique(dt$Seqs)))
    }

  } else {
    message("Task: extract specified-size sequences with input 'len' and 'step' under segment size cutoff.")
    dt$too_large <- (dt$end - dt$start + 1L) >= local_cutoff
    message("Total segments is ", nrow(dt), " and ", sum(dt$too_large), " of them with length >=", local_cutoff)
    message("Fraction: ", round(sum(dt$too_large) / nrow(dt), digits = 3))
    dt <- dt[, c("ID", "Seqs", "too_large")]
    dt <- dt[, list(Seqs = collapse_shift_seqs2(Seqs, too_large, len = len, step = step)),
             by = "ID"]

    if (return_dt) {
      return(dt)
    }

    all_seqs <- unique(dt$Seqs)
    keep <- nchar(all_seqs) >= len

    return(
      list(
        keep = sort(all_seqs[keep]),
        drop = sort(all_seqs[!keep])
      )
    )
  }
}

#' Get Copy Number Sequence Similarity or Distance Matrix
#'
#' @param x a coding copy number sequence (valid letters are A to X).
#' @param sub_mat default is `NULL`, use longest common substring method.
#' It can be a substitution matrix, each element indicates a score to plus.
#' See [sh_build_sub_matrix()].
#' @param block_size a block size to aggregrate, this is designed for big data, it means
#' results from adjacent sequences will be aggregrate by means to reduce the size of result
#' matrix.
#' @param dislike if `TRUE`, returns a dissimilarity matrix instead of a similarity matrix.
#' @param verbose if `TRUE`, print extra message, note it will slower the computation.
#' @param cores computer cores, default is `1`, note it is super fast already, set more
#' cores typically do not speed up the computation.
#' @inheritParams sh_build_sub_matrix
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
#' seqs2 <- sh_extract_seqs(x$dt, flexible_approach = TRUE)
#' seqs2
#'
#' mat <- sh_get_score_matrix(seqs$keep, x$mat, verbose = TRUE)
#' mat
#'
#' mat2 <- sh_get_score_matrix(seqs$keep, x$mat, dislike = TRUE)
#' identical(mat2, 120L - mat)
#'
#' mat_b <- sh_get_score_matrix(seqs$keep, x$mat, block_size = 2L)
#' ## block1 represents the first 2 sequences
#' ## block2 represents the 3rd, 4th sequences
#' ## ...
#' mat_b
#'
#' mat_c <- sh_get_score_matrix(seqs$keep)
#' mat_c
#' mat_d <- sh_get_score_matrix(seqs$keep, dislike = TRUE)
#' mat_d
#' \donttest{
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
#' }
#' @testexamples
#' expect_is(x, "list")
#' expect_is(seqs, "list")
#' expect_is(seqs2, "character")
#' expect_is(mat, "matrix")
#' expect_is(mat_b, "matrix")
#' expect_is(mat_c, "matrix")
#' expect_is(mat_d, "matrix")
#' expect_equal(mat2, 120L - mat)
#' if (require("doParallel")) {
#'   expect_equal(y1, y2)
#' }
sh_get_score_matrix <- function(x, sub_mat = NULL,
                                simple_version = FALSE,
                                block_size = NULL, dislike = FALSE,
                                cores = 1L, verbose = FALSE) {
  stopifnot(is.numeric(cores))

  if (is.null(sub_mat)) {
    message("Task: score paired strings with longest common substring method.")
    message("Final score = 2^(length - 1) for match otherwise 0.")
    message("  If dislike=TRUE, length of longest string - longest common substring is used.")

    y <- LCSMatrix(x, x, match = !dislike)
    y <- ifelse(y > 0, 2^(y - 1), 0)

    rownames(y) <- colnames(y) <- x

  } else {
    message("Task: score equal-size paired strings with substitution matrix.")
    if (anyNA(sub_mat)) {
      stop("Input substitution matrix cannot contain 'NA' values!")
    }

    if (simple_version) {
      map <- seq_len(6L)
      names(map) <- LETTERS[map]
    } else {
      map <- seq_len(24L)
      names(map) <- LETTERS[map]
    }
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
      y <- getScoreMatrix(m, sub_mat, block_size, !dislike, verbose)

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
        stop("Package 'foreach' is required to go through this parallel method.")
      }

      if (!"foreach" %in% .packages()) {
        attachNamespace("foreach")
      }

      if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("Package 'doParallel' is required to go through this parallel method.")
      }

      if (Sys.info()[["sysname"]] == "Windows") {
        cl <- parallel::makeCluster(cores)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl))
      } else {
        doParallel::registerDoParallel(cores = cores)
      }

      y <- foreach::foreach(
        i = seq_along(grp_list),
        .combine = "cbind",
        .packages = "sigminer.helper",
        # .export = c("m", "grp_list", "sub_mat", "verbose", "getScoreMatrixRect"),
        .export = c("getScoreMatrixRect"),
        .verbose = FALSE
      ) %dopar% {
        getScoreMatrixRect(m, m[grp_list[[i]], ], sub_mat, !dislike, verbose)
      }

      # for (i in seq_along(grp_list)) {
      #   y[, grp_list[[i]]] <- getScoreMatrixRect(m, m[grp_list[[i]], ], sub_mat, verbose)
      # }

      colnames(y) <- rownames(y) <- x
    }
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

#' Show Copy Number Sequence Shapes
#'
#' @inheritParams show_seq_logo
#' @inheritParams ggplot2::facet_wrap
#' @inheritParams sh_build_sub_matrix
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
#' x <- list(a = c("ABCDE", "AXFDP"), b = c("KKDFH", "GKDFM"))
#' p2 <- show_seq_shape(x)
#' p2
#' @testexamples
#' expect_is(p, "ggplot")
#' expect_is(p2, "ggplot")
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

#' Show Copy Number Sequence Logos
#' @inheritParams ggseqlogo2
#' @inheritParams sh_build_sub_matrix
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
