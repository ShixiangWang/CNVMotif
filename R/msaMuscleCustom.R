printHelp <- .get_fun("msa", "printHelp")
checkFunctionAvailable <- .get_fun("msa", "checkFunctionAvailable")
checkInputSeq <- .get_fun("msa", "checkInputSeq")
checkType <- .get_fun("msa", "checkType")
checkProfileScore <- .get_fun("msa", "checkProfileScore")
transformInputSeq <- .get_fun("msa", "transformInputSeq")
checkGapOpening2 <- .get_fun("msa", "checkGapOpening2")
checkGapExtension <- .get_fun("msa", "checkGapExtension")
checkMaxiters <- .get_fun("msa", "checkMaxiters")
checkVerbose <- .get_fun("msa", "checkVerbose")
checkIntegerParamsNew <- .get_fun("msa", "checkIntegerParamsNew")
checkNumericParamsNew <- .get_fun("msa", "checkNumericParamsNew")
checkNegativeParams <- .get_fun("msa", "checkNegativeParams")
checkSingleValParamsNew <- .get_fun("msa", "checkSingleValParamsNew")
checkPositiveParams <- .get_fun("msa", "checkPositiveParams")
checkIsValue <- .get_fun("msa", "checkIsValue")
checkInFile <- .get_fun("msa", "checkInFile")
checkIntervalParamsNew <- .get_fun("msa", "checkIntervalParamsNew")
checkLogicalParams <- .get_fun("msa", "checkLogicalParams")
convertAlnRows <- .get_fun("msa", "convertAlnRows")
printHelp <- .get_fun("msa", "printHelp")
printHelp <- .get_fun("msa", "printHelp")

msaMuscleCustom <- function(inputSeqs, cluster = "default", gapOpening = "default",
                            gapExtension = "default", maxiters = "default", substitutionMatrix = "default",
                            type = "default", order = c("aligned", "input"), verbose = FALSE,
                            help = FALSE, ...) {
  if (help) {
    printHelp("Muscle")
    return(invisible(NULL))
  }
  if (!checkFunctionAvailable("Muscle")) {
    stop("Muscle is not available via msa!")
  }
  params <- list(...)
  paramsCopy <- lapply(params, function(x) TRUE)
  params[["inputSeqIsFileFlag"]] <- checkInputSeq(inputSeqs)
  type <- checkType(type, inputSeqs, "msaMuscle")
  temporaryHelp <- checkProfileScore(type, params)
  params[["le"]] <- temporaryHelp[["le"]]
  params[["sp"]] <- temporaryHelp[["sp"]]
  params[["sv"]] <- temporaryHelp[["sv"]]
  params[["spn"]] <- temporaryHelp[["spn"]]
  inputSeqs <- transformInputSeq(inputSeqs)
  order <- match.arg(order)
  if (order == "input") {
    if (params[["inputSeqIsFileFlag"]]) {
      stop(
        "msaMuscle does not support order=\"input\" for reading\n",
        "sequences directly from a FASTA file."
      )
    } else if (is.null(names(inputSeqs)) || length(unique(names(inputSeqs))) !=
      length(inputSeqs)) {
      warning(
        "order=\"input\" requires input sequences to be named\n",
        "uniquely! Assigning default names 'Seq1'..'Seqn'\n",
        "to sequences."
      )
      names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))
    }
  }
  if (identical(cluster, "default") || is.null(cluster)) {
    cluster <- "upgma"
  }
  else {
    possibleValues <- c(
      "upgma", "upgmamax", "upgmamin",
      "upgmb", "neighborjoining"
    )
    if (length(cluster) != 1) {
      stop("The parameter cluster contains more than one value!")
    }
    if (!is.character(cluster)) {
      stop("The parameter cluster should be a string!")
    }
    cluster <- tolower(cluster)
    if (!(cluster %in% possibleValues)) {
      text <- ""
      text <- paste(possibleValues, collapse = ", ")
      stop(
        "The parameter cluster only can have the values: \n",
        text
      )
    }
  }
  if (is.null(substitutionMatrix) || identical(
    substitutionMatrix,
    "default"
  )) {
    substitutionMatrix <- NULL
  }
  if ((!is.null(substitutionMatrix) && !is.matrix(substitutionMatrix)) ||
    identical(mode(substitutionMatrix), "list")) {
    stop("The parameter substitutionMatrix should be a matrix!")
  }
  if (!is.null(substitutionMatrix)) {
    headerNames <- c(
      "A", "C", "D", "E", "F", "G", "H", "I",
      "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
      "W", "Y"
    )
    if (type == "protein") {
      reqNames <- headerNames
    } else if (type == "dna") {
      reqNames <- c("A", "C", "G", "T")
    } else {
      reqNames <- c("A", "C", "G", "U")
    }
    rowPerm <- match(reqNames, rownames(substitutionMatrix))
    if (any(is.na(rowPerm))) {
      stop("substitutionMatrix does not contain all necessary rows")
    }
    colPerm <- match(reqNames, colnames(substitutionMatrix))
    if (any(is.na(colPerm))) {
      stop("substitutionMatrix does not contain all necessary columns")
    }
    substitutionMatrix <- substitutionMatrix[rowPerm, colPerm]
    if (type == "rna") {
      reqNames <- c("A", "C", "G", "T")
    }
    auxMat <- matrix(0, length(headerNames), length(headerNames))
    rownames(auxMat) <- headerNames
    colnames(auxMat) <- headerNames
    auxMat[reqNames, reqNames] <- substitutionMatrix
    substitutionMatrix <- auxMat
    if (!isSymmetric(substitutionMatrix)) {
      stop("substitutionMatrix should be a symmetric matrix!")
    }
    if (any(is.na(substitutionMatrix)) || any(is.na(substitutionMatrix)) ||
      any(is.infinite(substitutionMatrix))) {
      stop("substitutionMatrix contains invalid values!")
    }
    params[["le"]] <- FALSE
    params[["sv"]] <- FALSE
    if (type == "protein") {
      params[["sp"]] <- TRUE
      params[["spn"]] <- FALSE
    }
    else {
      params[["sp"]] <- FALSE
      params[["spn"]] <- TRUE
    }
    paramsCopy[["le"]] <- NULL
    paramsCopy[["sv"]] <- NULL
    paramsCopy[["sp"]] <- NULL
    paramsCopy[["spn"]] <- NULL
  }
  if (params$le) {
    gapOpening <- checkGapOpening2(
      gapOpening, substitutionMatrix,
      2.9
    )
  }
  else if (params$sp) {
    gapOpening <- checkGapOpening2(
      gapOpening, substitutionMatrix,
      1439
    )
  }
  else if (params$sv) {
    gapOpening <- checkGapOpening2(
      gapOpening, substitutionMatrix,
      300
    )
  }
  else if (params$spn) {
    if (identical(type, "dna")) {
      gapOpening <- checkGapOpening2(
        gapOpening, substitutionMatrix,
        400
      )
    }
    if (identical(type, "rna")) {
      gapOpening <- checkGapOpening2(
        gapOpening, substitutionMatrix,
        420
      )
    }
    if (identical(type, "protein")) {
      stop(
        "If you use sequences of type \"protein\", \n",
        "you can't use the parameter \"spn\"!"
      )
    }
  }
  gapExtension <- checkGapExtension(
    gapExtension, type, substitutionMatrix,
    0, 0
  )
  maxiters <- checkMaxiters(maxiters, 16, "msaMuscle")
  verbose <- checkVerbose(FALSE, verbose)
  params[["anchorspacing"]] <- checkIntegerParamsNew(
    "anchorspacing",
    params
  )
  paramsCopy[["anchorspacing"]] <- NULL
  params[["center"]] <- checkNumericParamsNew("center", params)
  params[["center"]] <- checkNegativeParams("center", params)
  paramsCopy[["center"]] <- NULL
  posVal <- c("upgma", "upgmamax", "upgmamin", "upgmb", "neighborjoining")
  params[["cluster1"]] <- checkSingleValParamsNew(
    "cluster1",
    params, posVal
  )
  paramsCopy[["cluster1"]] <- NULL
  posVal <- c("upgma", "upgmb", "upgmamax", "upgmamin", "neighborjoining")
  params[["cluster2"]] <- checkSingleValParamsNew(
    "cluster2",
    params, posVal
  )
  paramsCopy[["cluster2"]] <- NULL
  params[["diagbreak"]] <- checkIntegerParamsNew(
    "diagbreak",
    params
  )
  params[["diagbreak"]] <- checkPositiveParams(
    "diagbreak",
    params
  )
  paramsCopy[["diagbreak"]] <- NULL
  params[["diaglength"]] <- checkIntegerParamsNew(
    "diaglength",
    params
  )
  params[["diaglength"]] <- checkPositiveParams(
    "diaglength",
    params
  )
  paramsCopy[["diaglength"]] <- NULL
  params[["diagmargin"]] <- checkIntegerParamsNew(
    "diagmargin",
    params
  )
  params[["diagmargin"]] <- checkPositiveParams(
    "diagmargin",
    params
  )
  paramsCopy[["diagmargin"]] <- NULL
  if (type == "protein") {
    posVal <- c("kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3")
  }
  else {
    posVal <- c(
      "kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3",
      "kmer4_6"
    )
  }
  if (!is.null(params[["distance1"]])) {
    params[["distance1"]] <- checkIsValue(
      "distance1", params,
      posVal
    )
  }
  paramsCopy[["distance1"]] <- NULL
  params[["distance2"]] <- checkSingleValParamsNew(
    "distance2",
    params, c("pctidkimura", "pctidlog")
  )
  paramsCopy[["distance2"]] <- NULL
  params[["hydro"]] <- checkIntegerParamsNew("hydro", params)
  paramsCopy[["hydro"]] <- NULL
  params[["hydrofactor"]] <- checkNumericParamsNew(
    "hydrofactor",
    params
  )
  paramsCopy[["hydrofactor"]] <- NULL
  if (!is.null(params[["in1"]])) {
    checkInFile("in1", params)
  }
  paramsCopy[["in1"]] <- NULL
  if (!is.null(params[["in2"]])) {
    checkInFile("in2", params)
  }
  paramsCopy[["in2"]] <- NULL
  params[["maxhours"]] <- checkNumericParamsNew(
    "maxhours",
    params
  )
  params[["maxhours"]] <- checkPositiveParams("maxhours", params)
  paramsCopy[["maxhours"]] <- NULL
  params[["maxtrees"]] <- checkIntegerParamsNew(
    "maxtrees",
    params
  )
  params[["maxtrees"]] <- checkPositiveParams("maxtrees", params)
  paramsCopy[["maxtrees"]] <- NULL
  params[["minbestcolscore"]] <- checkNumericParamsNew(
    "minbestcolscore",
    params
  )
  paramsCopy[["minbestcolscore"]] <- NULL
  params[["minsmoothscore"]] <- checkNumericParamsNew(
    "minsmoothscore",
    params
  )
  paramsCopy[["minsmoothscore"]] <- NULL
  posVal <- c("dp", "ps", "sp", "spf", "spm", "xp")
  params[["objscore"]] <- checkSingleValParamsNew(
    "objscore",
    params, posVal
  )
  paramsCopy[["objscore"]] <- NULL
  params[["refinewindow"]] <- checkIntegerParamsNew(
    "refinewindow",
    params
  )
  params[["refinewindow"]] <- checkPositiveParams(
    "refinewindow",
    params
  )
  paramsCopy[["refinewindow"]] <- NULL
  posVal <- c("pseudo", "midlongestspan", "minavgleafdist")
  params[["root1"]] <- checkSingleValParamsNew(
    "root1", params,
    posVal
  )
  paramsCopy[["root1"]] <- NULL
  posVal <- c("pseudo", "midlongestspan", "minavgleafdist")
  params[["root2"]] <- checkSingleValParamsNew(
    "root2", params,
    posVal
  )
  paramsCopy[["root2"]] <- NULL
  checkNumericParamsNew("smoothscoreceil", params)
  paramsCopy[["smoothscoreceil"]] <- NULL
  params[["smoothwindow"]] <- checkIntegerParamsNew(
    "smoothwindow",
    params
  )
  params[["smoothwindow"]] <- checkPositiveParams(
    "smoothwindow",
    params
  )
  if (!is.null(params[["smoothwindow"]]) && params[["smoothwindow"]] %% 2 ==
    0) {
    stop("The parameter smoothwindow must be odd!")
  }
  paramsCopy[["smoothwindow"]] <- NULL
  params[["SUEFF"]] <- checkNumericParamsNew("SUEFF", params)
  params[["SUEFF"]] <- checkIntervalParamsNew(
    "SUEFF", params,
    0, 1
  )
  paramsCopy[["SUEFF"]] <- NULL
  posVal <- c(
    "none", "henikoff", "henikoffpb", "gsc", "clustalw",
    "threeway"
  )
  params[["weight1"]] <- checkSingleValParamsNew(
    "weight1",
    params, posVal
  )
  paramsCopy[["weight1"]] <- NULL
  posVal <- c(
    "none", "henikoff", "henikoffpb", "gsc", "clustalw",
    "threeway"
  )
  params[["weight2"]] <- checkSingleValParamsNew(
    "weight2",
    params, posVal
  )
  paramsCopy[["weight2"]] <- NULL
  if (!is.null(params[["anchors"]])) {
    params[["anchors"]] <- checkLogicalParams(
      "anchors",
      params, TRUE
    )
  }
  paramsCopy[["anchors"]] <- NULL
  params[["brenner"]] <- checkLogicalParams(
    "brenner", params,
    FALSE
  )
  paramsCopy[["brenner"]] <- NULL
  if (!is.null(params[["core"]])) {
    params[["core"]] <- checkLogicalParams(
      "core", params,
      TRUE
    )
  }
  paramsCopy[["core"]] <- NULL
  params[["diags"]] <- checkLogicalParams(
    "diags", params,
    FALSE
  )
  paramsCopy[["diags"]] <- NULL
  params[["diags1"]] <- checkLogicalParams(
    "diags1", params,
    FALSE
  )
  paramsCopy[["diags1"]] <- NULL
  params[["diags2"]] <- checkLogicalParams(
    "diags2", params,
    FALSE
  )
  paramsCopy[["diags2"]] <- NULL
  params[["dimer"]] <- checkLogicalParams(
    "dimer", params,
    FALSE
  )
  paramsCopy[["dimer"]] <- NULL
  paramsCopy[["le"]] <- NULL
  params[["noanchors"]] <- checkLogicalParams(
    "noanchors",
    params, FALSE
  )
  if (!is.null(params[["anchors"]])) {
    if (params[["anchors"]] && params[["noanchors"]]) {
      stop("The parameters anchors and noanchors \n", "can't be positive at the same time!")
    }
    if (!params[["anchors"]] && !params[["noanchors"]]) {
      stop("The parameters anchors and noanchors \n", "can't be negative at the same time!")
    }
  }
  paramsCopy[["noanchors"]] <- NULL
  params[["nocore"]] <- checkLogicalParams(
    "nocore", params,
    FALSE
  )
  if (!is.null(params[["core"]])) {
    if (params[["core"]] && params[["nocore"]]) {
      stop("The parameters core and nocore \n", "can't be positive at the same time!")
    }
    if (!params[["core"]] && !params[["nocore"]]) {
      stop("The parameters core and nocore \n", "can't be negative at the same time!")
    }
  }
  paramsCopy[["nocore"]] <- NULL
  params[["profile"]] <- checkLogicalParams(
    "profile", params,
    FALSE
  )
  if (params[["profile"]]) {
    if (is.null(params[["in1"]]) || is.null(params[["in2"]])) {
      stop(
        "The parameter profile needs the following parameters: \n",
        "in1 and in2!"
      )
    }
  }
  paramsCopy[["profile"]] <- NULL
  params[["refine"]] <- checkLogicalParams(
    "refine", params,
    FALSE
  )
  paramsCopy[["refine"]] <- NULL
  params[["refinew"]] <- checkLogicalParams(
    "refinew", params,
    FALSE
  )
  paramsCopy[["refinew"]] <- NULL
  paramsCopy[["sp"]] <- NULL
  paramsCopy[["spn"]] <- NULL
  params[["spscore"]] <- checkLogicalParams(
    "spscore", params,
    FALSE
  )
  paramsCopy[["spscore"]] <- NULL
  paramsCopy[["sv"]] <- NULL
  params[["version"]] <- checkLogicalParams(
    "version", params,
    FALSE
  )
  paramsCopy[["version"]] <- NULL
  if (length(paramsCopy) != 0) {
    stop(
      "The following parameters are not known  \n", "(or have been specified",
      "more often than once):\n    ", paste(names(paramsCopy),
        collapse = ", ", sep = ""
      )
    )
  }
  inputSeqNames <- names(inputSeqs)
  names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))
  result <- .Call("RMuscle", inputSeqs, cluster, -abs(gapOpening),
    -abs(gapExtension), maxiters, substitutionMatrix, type,
    verbose, params,
    PACKAGE = "msa"
  )
  out <- convertAlnRows(result$msa, type)
  if (length(inputSeqNames) > 0) {
    if (order == "aligned") {
      perm <- match(names(out@unmasked), names(inputSeqs))
      names(out@unmasked) <- inputSeqNames[perm]
    }
    else {
      perm <- match(names(inputSeqs), names(out@unmasked))
      out@unmasked <- out@unmasked[perm]
      names(out@unmasked) <- inputSeqNames
    }
  }
  else {
    names(out@unmasked) <- NULL
  }
  standardParams <- list(
    gapOpening = gapOpening, gapExtension = gapExtension,
    maxiters = maxiters, verbose = verbose
  )
  out@params <- c(standardParams, params)
  out@call <- deparse(sys.call())
  out
}
