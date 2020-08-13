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
                                 seed = 1234L,
                                 ...) {
  stopifnot(is.matrix(x))

  message("Task: Estimate Optimal Number of Cluster for PAM Algorithm.")

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
    }
  } else if (method == "wss") {
    for (i in 1:k.max) {
      message("Generating ", i, " clusters...")
      clust <- FUNcluster(x, i, diss = TRUE, ...)
      v[i] <- .get_withinSS(diss, clust$cluster)
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
                      main = "Optimal number of clusters"
  )
  if (method == "silhouette") {
    p <- p + geom_vline(
      xintercept = which.max(v), linetype = 2,
      color = linecolor
    )
  }

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
