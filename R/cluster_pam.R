#' Estimate Optimal Number of Cluster for PAM Algorithm
#'
#' `cluster::clusGap()` cannot be used here for distance matrix, so
#' it is removed.
#'
#' @rdname cluster_pam
#' @inheritParams factoextra::fviz_nbclust
#' @inheritParams cluster::pam
#' @param x a dissimilarity matrix.
#' @param seed random seed.
#' @param clean_memory logical. If `TRUE`, the cluster result object will be removed
#' and the memory will be released by calling `gc()` to reduce the memory consumption.
#' @return a `ggplot` object.
#' @export
#' @examples
#' data("iris")
#' head(iris)
#' iris.scaled <- scale(iris[, -5])
#' iris.dist <- dist(iris.scaled) %>% as.matrix()
#' p <- cluster_pam_estimate(iris.dist)
#' p2 <- cluster_pam_estimate(iris.dist, method = "wss")
#'
#' cl <- cluster_pam(iris.dist, 3)
#' @testexamples
#' expect_is(p, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(cl, "pam")
cluster_pam_estimate <- function(x,
                                 method = c("silhouette", "wss"),
                                 k.max = 10, verbose = interactive(),
                                 barfill = "steelblue", barcolor = "steelblue", linecolor = "steelblue",
                                 FUNcluster = cluster::pam,
                                 seed = 1234L, clean_memory = FALSE,
                                 ...) {
  stopifnot(is.matrix(x))

  message("Task: Estimate Optimal Number of Cluster for PAM Algorithm.")
  on.exit(gc(verbose = FALSE))

  set.seed(seed)
  if (k.max < 2) {
    stop("k.max must bet > = 2")
  }
  method <- match.arg(method)

  diss <- stats::as.dist(x)
  v <- rep(0, k.max)
  if (method == "silhouette") {
    for (i in 2:k.max) {
      message("Generating ", i, " clusters...")
      clust <- FUNcluster(x, i, diss = TRUE, ...)
      v[i] <- .get_ave_sil_width(diss, clust$cluster)
      if (clean_memory) {
        rm(clust); gc(verbose = FALSE)
      }
    }
  } else if (method == "wss") {
    for (i in 1:k.max) {
      message("Generating ", i, " clusters...")
      clust <- FUNcluster(x, i, diss = TRUE, ...)
      v[i] <- .get_withinSS(diss, clust$cluster)
      if (clean_memory) {
        rm(clust); gc(verbose = FALSE)
      }
    }
  }
  message("Plotting...")
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
                      main = "Estimation of cluster number"
  )

  message("Done.")
  return(p)
}

# PAM Wrapper
#'
#' See [factoextra::fviz_cluster()] for visualizing result from `cluster_pam()`.
#' @inheritParams cluster::pam
#' @param x a dissimilarity matrix.
#' @param ... other parameters passing to [cluster::pam].
#' @return a PAM clustering result object.
#' @export
cluster_pam <- function(x, k, ...) {
  stopifnot(is.matrix(x))
  cluster::pam(x, k, diss = TRUE, ...)
}
