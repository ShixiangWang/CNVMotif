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

# Call un-exported functions from other packages
get_fun <- function(pkg, fun) {
  get(fun,
    envir = asNamespace(pkg),
    inherits = FALSE
  )
}

.get_ave_sil_width <- get_fun("factoextra", ".get_ave_sil_width")
.get_withinSS <- get_fun("factoextra", ".get_withinSS")

printHelp <- get_fun("msa", "printHelp")
checkFunctionAvailable <- get_fun("msa", "checkFunctionAvailable")
checkInputSeq <- get_fun("msa", "checkInputSeq")
checkType <- get_fun("msa", "checkType")
checkProfileScore <- get_fun("msa", "checkProfileScore")
transformInputSeq <- get_fun("msa", "transformInputSeq")
checkGapOpening2 <- get_fun("msa", "checkGapOpening2")
checkGapExtension <- get_fun("msa", "checkGapExtension")
checkMaxiters <- get_fun("msa", "checkMaxiters")
checkVerbose <- get_fun("msa", "checkVerbose")
checkIntegerParamsNew <- get_fun("msa", "checkIntegerParamsNew")
checkNumericParamsNew <- get_fun("msa", "checkNumericParamsNew")
checkNegativeParams <- get_fun("msa", "checkNegativeParams")
checkSingleValParamsNew <- get_fun("msa", "checkSingleValParamsNew")
checkPositiveParams <- get_fun("msa", "checkPositiveParams")
checkIsValue <- get_fun("msa", "checkIsValue")
checkInFile <- get_fun("msa", "checkInFile")
checkIntervalParamsNew <- get_fun("msa", "checkIntervalParamsNew")
checkLogicalParams <- get_fun("msa", "checkLogicalParams")
convertAlnRows <- get_fun("msa", "convertAlnRows")
checkGapOpening <- get_fun("msa", "checkGapOpening")
checkOutFile <- get_fun("msa", "checkOutFile")
interactiveCheck <- get_fun("msa", "interactiveCheck")
