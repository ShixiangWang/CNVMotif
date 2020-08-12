

#' Get Copy Number Sequence Similarity or Distance Matrix
#'
#' @param x a coding copy number sequence (valid letters are A to X).
#' @param sub_mat default is `NULL`, use longest common substring method.
#' It can be a substitution matrix, each element indicates a score to plus.
#' See [build_sub_matrix()].
#' @param block_size a block size to aggregrate, this is designed for big data, it means
#' results from adjacent sequences will be aggregrate by means to reduce the size of result
#' matrix.
#' @param dislike if `TRUE`, returns a dissimilarity matrix instead of a similarity matrix.
#' @param verbose if `TRUE`, print extra message, note it will slower the computation.
#' @param cores computer cores, default is `1`, note it is super fast already, set more
#' cores typically do not speed up the computation.
#' @inheritParams build_sub_matrix
#'
#' @return a score matrix.
#' @export
#'
#' @examples
#' load(system.file("extdata", "toy_segTab.RData",
#'   package = "CNVMotif", mustWork = TRUE
#' ))
#' x <- transform_seqs(segTabs)
#' x
#' seqs <- extract_seqs(x$dt)
#' seqs
#' seqs2 <- extract_seqs(x$dt, flexible_approach = TRUE)
#' seqs2
#'
#' mat <- get_score_matrix(seqs$keep, x$mat, verbose = TRUE)
#' mat
#'
#' mat2 <- get_score_matrix(seqs$keep, x$mat, dislike = TRUE)
#' identical(mat2, 120L - mat)
#'
#' mat_b <- get_score_matrix(seqs$keep, x$mat, block_size = 2L)
#' ## block1 represents the first 2 sequences
#' ## block2 represents the 3rd, 4th sequences
#' ## ...
#' mat_b
#'
#' mat_c <- get_score_matrix(seqs$keep)
#' mat_c
#' mat_d <- get_score_matrix(seqs$keep, dislike = TRUE)
#' mat_d
#' \donttest{
#' if (requireNamespace("doParallel")) {
#'   mock_seqs <- sapply(1:10000, function(x) {
#'     paste(sample(LETTERS[1:24], 5, replace = TRUE), collapse = "")
#'   })
#'
#'   system.time(
#'     y1 <- get_score_matrix(mock_seqs, x$mat, cores = 1)
#'   )
#'
#'   system.time(
#'     y2 <- get_score_matrix(mock_seqs, x$mat, cores = 2)
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
#' if (requireNamespace("doParallel")) {
#'   expect_equal(y1, y2)
#' }
get_score_matrix <- function(x, sub_mat = NULL,
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
        .packages = "CNVMotif",
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
