#' Extract Pasted Sequences from Each Chromosome
#'
#' See [get_score_matrix()] for examples.
#' The result sequences are unique and sorted.
#'
#' @param dt a `data.table` from [transform_seqs].
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
extract_seqs <- function(dt, len = 5L, step = 1L, local_cutoff = 1e7,
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
