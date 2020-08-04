#' Estimate Optimal Number of Cluster for PAM Algorithm
#' @inheritParams factoextra::fviz_nbclust
#' @inheritParams cluster::pam
#' @param x a dissimilarity matrix.
#' @return a `ggplot` object.
#' @export
cluster_pam_estimate <- function(x,
                                 method = c("silhouette", "wss", "gap_stat"),
                                 k.max = 10, nboot = 100, verbose = interactive(),
                                 barfill = "steelblue", barcolor = "steelblue", linecolor = "steelblue",
                                 print.summary = TRUE, FUNcluster = cluster::pam,
                                 seed = 1234L,
                                 ...) {
  stopifnot(is.matrix(x))
  set.seed(seed)
  if (k.max < 2) {
    stop("k.max must bet > = 2")
  }
  method <- match.arg(method)

  if (method %in% c("silhouette", "wss")) {
    diss <- stats::as.dist(x)
    v <- rep(0, k.max)
    if (method == "silhouette") {
      for (i in 2:k.max) {
        clust <- FUNcluster(x, i, diss = TRUE, ...)
        v[i] <- factoextra:::.get_ave_sil_width(diss, clust$cluster)
      }
    } else if (method == "wss") {
      for (i in 1:k.max) {
        clust <- FUNcluster(x, i, diss = TRUE, ...)
        v[i] <- factoextra:::.get_withinSS(diss, clust$cluster)
      }
    }
    df <- data.frame(
      clusters = as.factor(1:k.max), y = v,
      stringsAsFactors = TRUE
    )
    ylab <- "Total Within Sum of Square"
    if (method == "silhouette") {
      ylab <- "Average silhouette width"
    }
    p <- ggpubr::ggline(df,
      x = "clusters", y = "y", group = 1,
      color = linecolor, ylab = ylab, xlab = "Number of clusters k",
      main = "Optimal number of clusters"
    )
    if (method == "silhouette") {
      p <- p + geom_vline(
        xintercept = which.max(v), linetype = 2,
        color = linecolor
      )
    }
    return(p)
  } else if (method == "gap_stat") {
    extra_args <- list(...)
    gap_stat <- cluster::clusGap(x, FUNcluster,
      K.max = k.max,
      B = nboot, verbose = verbose,
      diss = TRUE,
      ...
    )
    if (!is.null(extra_args$maxSE)) {
      maxSE <- extra_args$maxSE
    } else {
      maxSE <- list(method = "firstSEmax", SE.factor = 1)
    }
    p <- factoextra::fviz_gap_stat(gap_stat, linecolor = linecolor, maxSE = maxSE)
    return(p)
  }
}

#' PAM Wrapper
#'
#' See [factoextra::fviz_cluster()] for visualization.
#' @inheritParams cluster::pam
#' @param x a dissimilarity matrix.
#' @return a PAM clustering result object.
#' @export
cluster_pam <- function(x, k, ...) {
  stopifnot(is.matrix(x))
  cluster::pam(x, k, diss = TRUE, ...)
}

