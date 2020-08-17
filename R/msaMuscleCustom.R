msaClustalWCustom <- function(inputSeqs, cluster = "default", gapOpening = "default",
                              gapExtension = "default", maxiters = "default",
                              substitutionMatrix = matrix(0),
                              type = "custom", order = c("input", "aligned"), verbose = FALSE,
                              help = FALSE, ...) {
  if (help) {
    printHelp("ClustalW")
    return(invisible(NULL))
  }
  if (!checkFunctionAvailable("ClustalW")) {
    stop("ClustalW is not available via msa!")
  }
  params <- list(...)
  paramsCopy <- lapply(params, function(x) TRUE)
  params[["inputSeqIsFileFlag"]] <- checkInputSeq(inputSeqs)
  # type <- checkType(type, inputSeqs, "msaClustalW")
  inputSeqs <- transformInputSeq(inputSeqs)
  order <- match.arg(order)
  params[["outorder"]] <- order
  if (is.null(cluster) || identical(cluster, "default")) {
    cluster <- "nj"
  }
  if (length(cluster) != 1) {
    stop(
      "The parameter cluster can only have one value. \n",
      "Possible values are \"nj\" or \"upgma\"!"
    )
  }
  cluster <- tolower(cluster)
  if (!(cluster %in% c("nj", "upgma"))) {
    stop("The parameter cluster can only have ", "the values \"nj\" or \"upgma\"!")
  }
  params[["substitutionMatrixIsDefaultFlag"]] <- FALSE
  params[["substitutionMatrixIsStringFlag"]] <- FALSE
  if (is.null(substitutionMatrix) || identical(
    substitutionMatrix,
    "default"
  )) {
    params[["substitutionMatrixIsDefaultFlag"]] <- TRUE
  }
  else if (is.character(substitutionMatrix) && !is.matrix(substitutionMatrix) &&
    grepl("\\.", substitutionMatrix, perl = TRUE)) {
    if (length(substitutionMatrix) != 1) {
      stop(
        "You are using more than one file for substitutionMatrix. \n",
        "It should only be a single character string!"
      )
    }
    if (!file.exists(substitutionMatrix)) {
      stop("The file for parameter substitutionMatrix does not exist!")
    }
    params[["substitutionMatrixIsStringFlag"]] <- TRUE
  }
  else if (is.character(substitutionMatrix) && !is.matrix(substitutionMatrix)) {
    if (type == "protein") {
      possibleValues <- c("blosum", "pam", "gonnet", "id")
      if (!(substitutionMatrix %in% possibleValues)) {
        text <- ""
        text <- paste(possibleValues, collapse = ", ")
        stop(
          "The parameter substitutionMatrix ", "only can have the values: \n",
          text
        )
      }
      params[["substitutionMatrixIsStringFlag"]] <- TRUE
    }
    else {
      possibleValues <- c("iub", "clustalw")
      if (!(substitutionMatrix %in% possibleValues)) {
        text <- ""
        text <- paste(possibleValues, collapse = ", ")
        stop(
          "The parameter substitutionMatrix ", "only can have the values: \n",
          text
        )
      }
      params[["substitutionMatrixIsStringFlag"]] <- FALSE
      params[["substitutionMatrixIsDefaultFlag"]] <- TRUE
      params[["pwdnamatrix"]] <- substitutionMatrix
      substitutionMatrix <- "default"
    }
  }
  else {
    reqNames <- c(
      "A", "R", "N", "D", "C", "Q", "E", "G",
      "H", "I", "L", "K", "M", "F", "P", "S", "T", "W",
      "Y", "V", "B", "Z", "X", "*"
    )
    if (type == "protein") {
      rowPerm <- match(reqNames, rownames(substitutionMatrix))
      if (any(is.na(rowPerm))) {
        stop("substitutionMatrix does not contain all necessary rows")
      }
      colPerm <- match(reqNames, colnames(substitutionMatrix))
      if (any(is.na(colPerm))) {
        stop("substitutionMatrix does not contain all necessary columns")
      }
      substitutionMatrix <- substitutionMatrix[
        rowPerm,
        colPerm
      ]
      if (!isSymmetric(substitutionMatrix)) {
        stop("substitutionMatrix should be a symmetric matrix!")
      }
    }
    else if (type %in% c("dna", "rna")) {
      reqNuc <- if (type == "dna") {
        c("A", "G", "C", "T")
      } else {
        c("A", "G", "C", "U")
      }
      if (any(is.na(match(reqNuc, rownames(substitutionMatrix))))) {
        stop("substitutionMatrix does not contain all necessary rows")
      }
      if (any(is.na(match(reqNuc, colnames(substitutionMatrix))))) {
        stop("substitutionMatrix does not contain all necessary columns")
      }
      rowSel <- which(rownames(substitutionMatrix) %in%
        reqNames)
      colSel <- which(colnames(substitutionMatrix) %in%
        reqNames)
      substitutionMatrix <- substitutionMatrix[
        rowSel,
        colSel
      ]
      fakeAAmat <- matrix(0, length(reqNames), length(reqNames))
      rownames(fakeAAmat) <- reqNames
      colnames(fakeAAmat) <- reqNames
      fakeAAmat[rownames(substitutionMatrix), colnames(substitutionMatrix)] <- substitutionMatrix
      substitutionMatrix <- fakeAAmat
      params[["dnamatrix"]] <- NULL
    } else {
      ## TO support custom letters
      type <- "protein"
      reqNames <- c(setdiff(LETTERS, c("J")), "*")
      if (!all(reqNames %in% rownames(substitutionMatrix))) {
        substitutionMatrix2 <- matrix(0,
          nrow = length(reqNames), ncol = length(reqNames),
          dimnames = list(reqNames, reqNames)
        )
        diag(substitutionMatrix2) <- 1
        substitutionMatrix2["*", "*"] <- 0
        index <- intersect(rownames(substitutionMatrix), rownames(substitutionMatrix2))
        substitutionMatrix2[index, index] <- substitutionMatrix[index, index]
        substitutionMatrix <- substitutionMatrix2
      }
    }
  }
  gapOpening <- checkGapOpening(gapOpening, type, substitutionMatrix,
    defaultDNAValue = 15, defaultAAValue = 10
  )
  gapExtension <- checkGapExtension(gapExtension, type, substitutionMatrix,
    defaultDNAValue = 6.66, defaultAAValue = 0.2
  )
  maxiters <- checkMaxiters(maxiters, 3, "msaClustalW")
  verbose <- checkVerbose(FALSE, verbose)
  params[["options"]] <- checkLogicalParams(
    "options", params,
    FALSE
  )
  paramsCopy[["options"]] <- NULL
  params[["check"]] <- checkLogicalParams(
    "check", params,
    FALSE
  )
  paramsCopy[["check"]] <- NULL
  params[["fullhelp"]] <- checkLogicalParams(
    "fullhelp", params,
    FALSE
  )
  paramsCopy[["fullhelp"]] <- NULL
  params[["align"]] <- checkLogicalParams(
    "align", params,
    FALSE
  )
  paramsCopy[["align"]] <- NULL
  params[["pim"]] <- checkLogicalParams("pim", params, FALSE)
  paramsCopy[["pim"]] <- NULL
  params[["convert"]] <- checkLogicalParams(
    "convert", params,
    FALSE
  )
  paramsCopy[["convert"]] <- NULL
  params[["quicktree"]] <- checkLogicalParams(
    "quicktree",
    params, FALSE
  )
  paramsCopy[["quicktree"]] <- NULL
  params[["negative"]] <- checkLogicalParams(
    "negative", params,
    FALSE
  )
  paramsCopy[["negative"]] <- NULL
  posVal <- "clustal"
  if (!is.null(params[["output"]]) && !identical(
    params[["output"]],
    posVal
  )) {
    stop(
      "Until now, the only value for parameter \n", "output is \"clustal\", which is default. \n",
      "A more sophisticated implementation should \n",
      "be available in higher versions of the package."
    )
  }
  params[["output"]] <- checkSingleValParamsNew(
    "output", params,
    posVal
  )
  paramsCopy[["output"]] <- NULL
  posVal <- c("lower", "upper")
  params[["case"]] <- checkSingleValParamsNew(
    "case", params,
    posVal
  )
  paramsCopy[["case"]] <- NULL
  if (is.null(params[["seqnos"]])) {
    params[["seqnosFlag"]] <- TRUE
  }
  else {
    params[["seqnosFlag"]] <- FALSE
    if (length(params[["seqnos"]]) != 1) {
      stop(
        "The parameter seqnos should be a single string! \n",
        "Possible values are \"on\", or \"off\"!"
      )
    }
    if (!is.character(params[["seqnos"]])) {
      stop(
        "The parameter <seqnos> should be a string! \n",
        "Possible values are \"on\", or \"off\"!"
      )
    }
    posVal <- c("on", "off")
    params[["seqnos"]] <- checkIsValue(
      "seqnos", params,
      posVal
    )
  }
  paramsCopy[["seqnos"]] <- NULL
  posVal <- c("off", "on")
  params[["seqno_range"]] <- checkSingleValParamsNew(
    "seqno_range",
    params, posVal
  )
  paramsCopy[["seqno_range"]] <- NULL
  if (!is.null(params[["range"]])) {
    if (length(params[["range"]]) != 2) {
      stop(
        "The parameter range needs a vector of length 2! \n",
        "Both values should be positive integers!"
      )
    }
    if (!is.vector(params[["range"]])) {
      stop("The parameter range should be a vector ", "with 2 positive integers!")
    }
    if (any(is.na(params[["range"]])) || any(is.nan(params[["range"]]))) {
      stop(
        "The parameter range should consist of 2 positive ",
        "integers, \n", "not with NAs or NaNs!"
      )
    }
    if (!is.integer(params[["range"]])) {
      if (params[["range"]][[1]] - round(params[["range"]][[1]]) !=
        0 | params[["range"]][[2]] - round(params[["range"]][[2]]) !=
        0) {
        stop(
          "The parameter range should consist of integers, \n",
          "not numeric values!"
        )
      }
      if (params[["range"]][[1]] <= .Machine$integer.max &
        params[["range"]][[2]] <= .Machine$integer.max) {
        params[["range"]][[1]] <- as.integer(params[["range"]][[1]])
        params[["range"]][[2]] <- as.integer(params[["range"]][[2]])
      }
      else {
        stop("The values in parameter range ", " are bigger than integer!")
      }
    }
    if (params[["range"]][[1]] < 0 | params[["range"]][[2]] <
      0) {
      stop("The parameter range needs positive integer values!")
    }
  }
  paramsCopy[["range"]] <- NULL
  if (!is.null(params[["stats"]])) {
    tempList <- checkOutFile("stats", params)
    if (tempList$existingFile) {
      interactiveCheck("stats", params)
    }
    params[["stats"]] <- tempList$param
  }
  paramsCopy[["stats"]] <- NULL
  params[["ktuple"]] <- checkIntegerParamsNew("ktuple", params)
  if (!is.null(params[["ktuple"]])) {
    if (type == "protein") {
      if (params[["ktuple"]] > 2) {
        stop("If you are using proteins, ktuple should be <=2!")
      }
    }
  }
  paramsCopy[["ktuple"]] <- NULL
  params[["topdiags"]] <- checkIntegerParamsNew(
    "topdiags",
    params
  )
  paramsCopy[["topdiags"]] <- NULL
  params[["window"]] <- checkIntegerParamsNew("window", params)
  paramsCopy[["window"]] <- NULL
  params[["pairgap"]] <- checkIntegerParamsNew("pairgap", params)
  paramsCopy[["pairgap"]] <- NULL
  posVal <- c("percent", "absolute")
  params[["score"]] <- checkSingleValParamsNew(
    "score", params,
    posVal
  )
  paramsCopy[["score"]] <- NULL
  if (!is.null(params[["pwmatrix"]]) && grepl("\\.", params[["pwmatrix"]],
    perl = TRUE
  )) {
    checkInFile("pwmatrix", params)
  }
  else {
    posVal <- c("blosum", "pam", "gonnet", "id")
    params[["pwmatrix"]] <- checkSingleValParamsNew(
      "pwmatrix",
      params, posVal
    )
  }
  paramsCopy[["pwmatrix"]] <- NULL
  if (!is.null(params[["pwdnamatrix"]]) && grepl("\\.", params[["pwdnamatrix"]],
    perl = TRUE
  )) {
    checkInFile("pwdnamatrix", params)
  }
  else {
    posVal <- c("iub", "clustalw")
    params[["pwdnamatrix"]] <- checkSingleValParamsNew(
      "pwdnamatrix",
      params, posVal
    )
  }
  paramsCopy[["pwdnamatrix"]] <- NULL
  params[["pwgapopen"]] <- checkNumericParamsNew(
    "pwgapopen",
    params
  )
  if (is.numeric(params[["pwgapopen"]])) {
    params[["pwgapopen"]] <- abs(params[["pwgapopen"]])
  }
  paramsCopy[["pwgapopen"]] <- NULL
  params[["pwgapext"]] <- checkNumericParamsNew(
    "pwgapext",
    params
  )
  if (is.numeric(params[["pwgapext"]])) {
    params[["pwgapext"]] <- abs(params[["pwgapext"]])
  }
  paramsCopy[["pwgapext"]] <- NULL
  if (!is.null(params[["usetree"]])) {
    checkInFile("usetree", params)
  }
  paramsCopy[["usetree"]] <- NULL
  if (!is.null(params[["dnamatrix"]]) && grepl("\\.", params[["dnamatrix"]],
    perl = TRUE
  )) {
    checkInFile("dnamatrix", params)
  }
  else if (is.null(params[["pwdnamatrix"]])) {
    posVal <- c("iub", "clustalw")
    params[["pwdnamatrix"]] <- checkSingleValParamsNew(
      "dnamatrix",
      params, posVal
    )
  }
  paramsCopy[["dnamatrix"]] <- NULL
  params[["endgaps"]] <- checkLogicalParams(
    "endgaps", params,
    FALSE
  )
  paramsCopy[["endgaps"]] <- NULL
  params[["gapdist"]] <- checkIntegerParamsNew("gapdist", params)
  paramsCopy[["gapdist"]] <- NULL
  params[["nopgap"]] <- checkLogicalParams(
    "nopgap", params,
    FALSE
  )
  paramsCopy[["nopgap"]] <- NULL
  params[["nohgap"]] <- checkLogicalParams(
    "nohgap", params,
    FALSE
  )
  paramsCopy[["nohgap"]] <- NULL
  if (!is.null(params[["novgap"]])) {
    params[["novgap"]] <- checkLogicalParams(
      "novgap", params,
      TRUE
    )
  }
  paramsCopy[["novgap"]] <- NULL
  if (!is.null(params[["hgapresidues"]])) {
    if (!is.character(params[["hgapresidues"]])) {
      stop("The parameter hgapresidues should be a string!")
    }
  }
  paramsCopy[["hgapresidues"]] <- NULL
  params[["maxdiv"]] <- checkIntegerParamsNew("maxdiv", params)
  paramsCopy[["maxdiv"]] <- NULL
  params[["transweight"]] <- checkNumericParamsNew(
    "transweight",
    params
  )
  paramsCopy[["transweight"]] <- NULL
  posVal <- c("tree", "alignment", "none")
  params[["iteration"]] <- checkSingleValParamsNew(
    "iteration",
    params, posVal
  )
  paramsCopy[["iteration"]] <- NULL
  params[["noweights"]] <- checkLogicalParams(
    "noweights",
    params, FALSE
  )
  paramsCopy[["noweights"]] <- NULL
  params[["profile"]] <- checkLogicalParams(
    "profile", params,
    FALSE
  )
  paramsCopy[["profile"]] <- NULL
  if (!is.null(params[["profile1"]])) {
    checkInFile("profile1", params)
  }
  paramsCopy[["profile1"]] <- NULL
  if (!is.null(params[["profile2"]])) {
    checkInFile("profile2", params)
  }
  paramsCopy[["profile2"]] <- NULL
  if (!is.null(params[["usetree1"]])) {
    checkInFile("usetree1", params)
  }
  paramsCopy[["usetree1"]] <- NULL
  if (!is.null(params[["usetree2"]])) {
    checkInFile("usetree2", params)
  }
  paramsCopy[["usetree2"]] <- NULL
  params[["sequences"]] <- checkLogicalParams(
    "sequences",
    params, FALSE
  )
  paramsCopy[["sequences"]] <- NULL
  params[["nosecstr1"]] <- checkLogicalParams(
    "nosecstr1",
    params, FALSE
  )
  paramsCopy[["nosecstr1"]] <- NULL
  params[["nosecstr2"]] <- checkLogicalParams(
    "nosecstr2",
    params, FALSE
  )
  paramsCopy[["nosecstr2"]] <- NULL
  posVal <- c("structure", "mask", "both", "none")
  params[["secstrout"]] <- checkSingleValParamsNew(
    "secstrout",
    params, posVal
  )
  paramsCopy[["secstrout"]] <- NULL
  params[["helixgap"]] <- checkIntegerParamsNew(
    "helixgap",
    params
  )
  paramsCopy[["helixgap"]] <- NULL
  params[["strandgap"]] <- checkIntegerParamsNew(
    "strandgap",
    params
  )
  paramsCopy[["strandgap"]] <- NULL
  params[["loopgap"]] <- checkIntegerParamsNew("loopgap", params)
  paramsCopy[["loopgap"]] <- NULL
  params[["terminalgap"]] <- checkIntegerParamsNew(
    "terminalgap",
    params
  )
  paramsCopy[["terminalgap"]] <- NULL
  params[["helixendin"]] <- checkIntegerParamsNew(
    "helixendin",
    params
  )
  paramsCopy[["helixendin"]] <- NULL
  params[["helixendout"]] <- checkIntegerParamsNew(
    "helixendout",
    params
  )
  paramsCopy[["helixendout"]] <- NULL
  params[["strandendin"]] <- checkIntegerParamsNew(
    "strandendin",
    params
  )
  paramsCopy[["strandendin"]] <- NULL
  params[["strandendout"]] <- checkIntegerParamsNew(
    "strandendout",
    params
  )
  paramsCopy[["strandendout"]] <- NULL
  posVal <- c("nj", "phylip", "dist", "nexus")
  params[["outputtree"]] <- checkSingleValParamsNew(
    "outputtree",
    params, posVal
  )
  paramsCopy[["outputtree"]] <- NULL
  params[["seed"]] <- checkIntegerParamsNew("seed", params)
  paramsCopy[["seed"]] <- NULL
  params[["kimura"]] <- checkLogicalParams(
    "kimura", params,
    FALSE
  )
  paramsCopy[["kimura"]] <- NULL
  params[["tossgaps"]] <- checkLogicalParams(
    "tossgaps", params,
    FALSE
  )
  paramsCopy[["tossgaps"]] <- NULL
  posVal <- c("node", "branch")
  params[["bootlabels"]] <- checkSingleValParamsNew(
    "bootlabels",
    params, posVal
  )
  paramsCopy[["bootlabels"]] <- NULL
  if (length(paramsCopy) != 0) {
    stop(
      "The following parameters are not known \n", "(or have been specified",
      "more often than once):\n    ", paste(names(paramsCopy),
        collapse = ", ", sep = ""
      )
    )
  }
  inputSeqNames <- names(inputSeqs)
  names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))
  result <- .Call("RClustalW", inputSeqs, cluster, abs(gapOpening),
    abs(gapExtension), maxiters, substitutionMatrix, type,
    verbose, params,
    PACKAGE = "msa"
  )
  out <- convertAlnRows(result$msa, type)
  if (length(inputSeqNames) > 0) {
    perm <- match(names(out@unmasked), names(inputSeqs))
    names(out@unmasked) <- inputSeqNames[perm]
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


# msaMuscleCustom <- function(inputSeqs, cluster = "default", gapOpening = "default",
#                             gapExtension = "default", maxiters = "default", substitutionMatrix = matrix(0),
#                             type = "custom", order = c("input", "aligned"), verbose = FALSE,
#                             help = FALSE, ...) {
#   if (help) {
#     printHelp("Muscle")
#     return(invisible(NULL))
#   }
#   if (!checkFunctionAvailable("Muscle")) {
#     stop("Muscle is not available via msa!")
#   }
#   params <- list(...)
#   paramsCopy <- lapply(params, function(x) TRUE)
#   params[["inputSeqIsFileFlag"]] <- checkInputSeq(inputSeqs)
#   # type <- checkType(type, inputSeqs, "msaMuscle")
#   temporaryHelp <- checkProfileScore(type, params)
#   params[["le"]] <- temporaryHelp[["le"]]
#   params[["sp"]] <- temporaryHelp[["sp"]]
#   params[["sv"]] <- temporaryHelp[["sv"]]
#   params[["spn"]] <- temporaryHelp[["spn"]]
#   inputSeqs <- transformInputSeq(inputSeqs)
#   order <- match.arg(order)
#   if (order == "input") {
#     if (params[["inputSeqIsFileFlag"]]) {
#       stop(
#         "msaMuscle does not support order=\"input\" for reading\n",
#         "sequences directly from a FASTA file."
#       )
#     } else if (is.null(names(inputSeqs)) || length(unique(names(inputSeqs))) !=
#       length(inputSeqs)) {
#       warning(
#         "order=\"input\" requires input sequences to be named\n",
#         "uniquely! Assigning default names 'Seq1'..'Seqn'\n",
#         "to sequences."
#       )
#       names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))
#     }
#   }
#   if (identical(cluster, "default") || is.null(cluster)) {
#     cluster <- "upgma"
#   }
#   else {
#     possibleValues <- c(
#       "upgma", "upgmamax", "upgmamin",
#       "upgmb", "neighborjoining"
#     )
#     if (length(cluster) != 1) {
#       stop("The parameter cluster contains more than one value!")
#     }
#     if (!is.character(cluster)) {
#       stop("The parameter cluster should be a string!")
#     }
#     cluster <- tolower(cluster)
#     if (!(cluster %in% possibleValues)) {
#       text <- ""
#       text <- paste(possibleValues, collapse = ", ")
#       stop(
#         "The parameter cluster only can have the values: \n",
#         text
#       )
#     }
#   }
#   if (is.null(substitutionMatrix) || identical(
#     substitutionMatrix,
#     "default"
#   )) {
#     substitutionMatrix <- NULL
#   }
#   if ((!is.null(substitutionMatrix) && !is.matrix(substitutionMatrix)) ||
#     identical(mode(substitutionMatrix), "list")) {
#     stop("The parameter substitutionMatrix should be a matrix!")
#   }
#   if (!is.null(substitutionMatrix)) {
#     headerNames <- c(
#       "A", "C", "D", "E", "F", "G", "H", "I",
#       "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
#       "W", "Y"
#     )
#     if (type == "protein") {
#       reqNames <- headerNames
#     } else if (type == "dna") {
#       reqNames <- c("A", "C", "G", "T")
#     } else if (type == "rna") {
#       reqNames <- c("A", "C", "G", "U")
#     }
#
#     if (!type %in% c("dna", "rna", "protein")) {
#       ## TO support custom letters
#       type <- "protein"
#       reqNames <- setdiff(LETTERS, c("J", "O", "U"))
#       if (!all(reqNames %in% rownames(substitutionMatrix))) {
#         substitutionMatrix2 <- matrix(0,
#                                       nrow = length(reqNames), ncol = length(reqNames),
#                                       dimnames = list(reqNames, reqNames))
#         diag(substitutionMatrix2) <- 1
#         index <- intersect(rownames(substitutionMatrix), rownames(substitutionMatrix2))
#         substitutionMatrix2[index, index] <- substitutionMatrix[index, index]
#         substitutionMatrix <- substitutionMatrix2
#       }
#     } else {
#       rowPerm <- match(reqNames, rownames(substitutionMatrix))
#       if (any(is.na(rowPerm))) {
#         stop("substitutionMatrix does not contain all necessary rows")
#       }
#       colPerm <- match(reqNames, colnames(substitutionMatrix))
#       if (any(is.na(colPerm))) {
#         stop("substitutionMatrix does not contain all necessary columns")
#       }
#       substitutionMatrix <- substitutionMatrix[rowPerm, colPerm]
#       if (type == "rna") {
#         reqNames <- c("A", "C", "G", "T")
#       }
#       auxMat <- matrix(0, length(headerNames), length(headerNames))
#       rownames(auxMat) <- headerNames
#       colnames(auxMat) <- headerNames
#       auxMat[reqNames, reqNames] <- substitutionMatrix
#       substitutionMatrix <- auxMat
#     }
#
#     if (!isSymmetric(substitutionMatrix)) {
#       stop("substitutionMatrix should be a symmetric matrix!")
#     }
#     if (any(is.na(substitutionMatrix)) || any(is.na(substitutionMatrix)) ||
#       any(is.infinite(substitutionMatrix))) {
#       stop("substitutionMatrix contains invalid values!")
#     }
#     params[["le"]] <- FALSE
#     params[["sv"]] <- FALSE
#     if (type == "protein") {
#       params[["sp"]] <- TRUE
#       params[["spn"]] <- FALSE
#     }
#     else {
#       params[["sp"]] <- FALSE
#       params[["spn"]] <- TRUE
#     }
#     paramsCopy[["le"]] <- NULL
#     paramsCopy[["sv"]] <- NULL
#     paramsCopy[["sp"]] <- NULL
#     paramsCopy[["spn"]] <- NULL
#   }
#   if (params$le) {
#     gapOpening <- checkGapOpening2(
#       gapOpening, substitutionMatrix,
#       2.9
#     )
#   }
#   else if (params$sp) {
#     gapOpening <- checkGapOpening2(
#       gapOpening, substitutionMatrix,
#       1439
#     )
#   }
#   else if (params$sv) {
#     gapOpening <- checkGapOpening2(
#       gapOpening, substitutionMatrix,
#       300
#     )
#   }
#   else if (params$spn) {
#     if (identical(type, "dna")) {
#       gapOpening <- checkGapOpening2(
#         gapOpening, substitutionMatrix,
#         400
#       )
#     }
#     if (identical(type, "rna")) {
#       gapOpening <- checkGapOpening2(
#         gapOpening, substitutionMatrix,
#         420
#       )
#     }
#     if (identical(type, "protein")) {
#       stop(
#         "If you use sequences of type \"protein\", \n",
#         "you can't use the parameter \"spn\"!"
#       )
#     }
#   }
#   gapExtension <- checkGapExtension(
#     gapExtension, type, substitutionMatrix,
#     0, 0
#   )
#   maxiters <- checkMaxiters(maxiters, 16, "msaMuscle")
#   verbose <- checkVerbose(FALSE, verbose)
#   params[["anchorspacing"]] <- checkIntegerParamsNew(
#     "anchorspacing",
#     params
#   )
#   paramsCopy[["anchorspacing"]] <- NULL
#   params[["center"]] <- checkNumericParamsNew("center", params)
#   params[["center"]] <- checkNegativeParams("center", params)
#   paramsCopy[["center"]] <- NULL
#   posVal <- c("upgma", "upgmamax", "upgmamin", "upgmb", "neighborjoining")
#   params[["cluster1"]] <- checkSingleValParamsNew(
#     "cluster1",
#     params, posVal
#   )
#   paramsCopy[["cluster1"]] <- NULL
#   posVal <- c("upgma", "upgmb", "upgmamax", "upgmamin", "neighborjoining")
#   params[["cluster2"]] <- checkSingleValParamsNew(
#     "cluster2",
#     params, posVal
#   )
#   paramsCopy[["cluster2"]] <- NULL
#   params[["diagbreak"]] <- checkIntegerParamsNew(
#     "diagbreak",
#     params
#   )
#   params[["diagbreak"]] <- checkPositiveParams(
#     "diagbreak",
#     params
#   )
#   paramsCopy[["diagbreak"]] <- NULL
#   params[["diaglength"]] <- checkIntegerParamsNew(
#     "diaglength",
#     params
#   )
#   params[["diaglength"]] <- checkPositiveParams(
#     "diaglength",
#     params
#   )
#   paramsCopy[["diaglength"]] <- NULL
#   params[["diagmargin"]] <- checkIntegerParamsNew(
#     "diagmargin",
#     params
#   )
#   params[["diagmargin"]] <- checkPositiveParams(
#     "diagmargin",
#     params
#   )
#   paramsCopy[["diagmargin"]] <- NULL
#   if (type == "protein") {
#     posVal <- c("kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3")
#   }
#   else {
#     posVal <- c(
#       "kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3",
#       "kmer4_6"
#     )
#   }
#   if (!is.null(params[["distance1"]])) {
#     params[["distance1"]] <- checkIsValue(
#       "distance1", params,
#       posVal
#     )
#   }
#   paramsCopy[["distance1"]] <- NULL
#   params[["distance2"]] <- checkSingleValParamsNew(
#     "distance2",
#     params, c("pctidkimura", "pctidlog")
#   )
#   paramsCopy[["distance2"]] <- NULL
#   params[["hydro"]] <- checkIntegerParamsNew("hydro", params)
#   paramsCopy[["hydro"]] <- NULL
#   params[["hydrofactor"]] <- checkNumericParamsNew(
#     "hydrofactor",
#     params
#   )
#   paramsCopy[["hydrofactor"]] <- NULL
#   if (!is.null(params[["in1"]])) {
#     checkInFile("in1", params)
#   }
#   paramsCopy[["in1"]] <- NULL
#   if (!is.null(params[["in2"]])) {
#     checkInFile("in2", params)
#   }
#   paramsCopy[["in2"]] <- NULL
#   params[["maxhours"]] <- checkNumericParamsNew(
#     "maxhours",
#     params
#   )
#   params[["maxhours"]] <- checkPositiveParams("maxhours", params)
#   paramsCopy[["maxhours"]] <- NULL
#   params[["maxtrees"]] <- checkIntegerParamsNew(
#     "maxtrees",
#     params
#   )
#   params[["maxtrees"]] <- checkPositiveParams("maxtrees", params)
#   paramsCopy[["maxtrees"]] <- NULL
#   params[["minbestcolscore"]] <- checkNumericParamsNew(
#     "minbestcolscore",
#     params
#   )
#   paramsCopy[["minbestcolscore"]] <- NULL
#   params[["minsmoothscore"]] <- checkNumericParamsNew(
#     "minsmoothscore",
#     params
#   )
#   paramsCopy[["minsmoothscore"]] <- NULL
#   posVal <- c("dp", "ps", "sp", "spf", "spm", "xp")
#   params[["objscore"]] <- checkSingleValParamsNew(
#     "objscore",
#     params, posVal
#   )
#   paramsCopy[["objscore"]] <- NULL
#   params[["refinewindow"]] <- checkIntegerParamsNew(
#     "refinewindow",
#     params
#   )
#   params[["refinewindow"]] <- checkPositiveParams(
#     "refinewindow",
#     params
#   )
#   paramsCopy[["refinewindow"]] <- NULL
#   posVal <- c("pseudo", "midlongestspan", "minavgleafdist")
#   params[["root1"]] <- checkSingleValParamsNew(
#     "root1", params,
#     posVal
#   )
#   paramsCopy[["root1"]] <- NULL
#   posVal <- c("pseudo", "midlongestspan", "minavgleafdist")
#   params[["root2"]] <- checkSingleValParamsNew(
#     "root2", params,
#     posVal
#   )
#   paramsCopy[["root2"]] <- NULL
#   checkNumericParamsNew("smoothscoreceil", params)
#   paramsCopy[["smoothscoreceil"]] <- NULL
#   params[["smoothwindow"]] <- checkIntegerParamsNew(
#     "smoothwindow",
#     params
#   )
#   params[["smoothwindow"]] <- checkPositiveParams(
#     "smoothwindow",
#     params
#   )
#   if (!is.null(params[["smoothwindow"]]) && params[["smoothwindow"]] %% 2 ==
#     0) {
#     stop("The parameter smoothwindow must be odd!")
#   }
#   paramsCopy[["smoothwindow"]] <- NULL
#   params[["SUEFF"]] <- checkNumericParamsNew("SUEFF", params)
#   params[["SUEFF"]] <- checkIntervalParamsNew(
#     "SUEFF", params,
#     0, 1
#   )
#   paramsCopy[["SUEFF"]] <- NULL
#   posVal <- c(
#     "none", "henikoff", "henikoffpb", "gsc", "clustalw",
#     "threeway"
#   )
#   params[["weight1"]] <- checkSingleValParamsNew(
#     "weight1",
#     params, posVal
#   )
#   paramsCopy[["weight1"]] <- NULL
#   posVal <- c(
#     "none", "henikoff", "henikoffpb", "gsc", "clustalw",
#     "threeway"
#   )
#   params[["weight2"]] <- checkSingleValParamsNew(
#     "weight2",
#     params, posVal
#   )
#   paramsCopy[["weight2"]] <- NULL
#   if (!is.null(params[["anchors"]])) {
#     params[["anchors"]] <- checkLogicalParams(
#       "anchors",
#       params, TRUE
#     )
#   }
#   paramsCopy[["anchors"]] <- NULL
#   params[["brenner"]] <- checkLogicalParams(
#     "brenner", params,
#     FALSE
#   )
#   paramsCopy[["brenner"]] <- NULL
#   if (!is.null(params[["core"]])) {
#     params[["core"]] <- checkLogicalParams(
#       "core", params,
#       TRUE
#     )
#   }
#   paramsCopy[["core"]] <- NULL
#   params[["diags"]] <- checkLogicalParams(
#     "diags", params,
#     FALSE
#   )
#   paramsCopy[["diags"]] <- NULL
#   params[["diags1"]] <- checkLogicalParams(
#     "diags1", params,
#     FALSE
#   )
#   paramsCopy[["diags1"]] <- NULL
#   params[["diags2"]] <- checkLogicalParams(
#     "diags2", params,
#     FALSE
#   )
#   paramsCopy[["diags2"]] <- NULL
#   params[["dimer"]] <- checkLogicalParams(
#     "dimer", params,
#     FALSE
#   )
#   paramsCopy[["dimer"]] <- NULL
#   paramsCopy[["le"]] <- NULL
#   params[["noanchors"]] <- checkLogicalParams(
#     "noanchors",
#     params, FALSE
#   )
#   if (!is.null(params[["anchors"]])) {
#     if (params[["anchors"]] && params[["noanchors"]]) {
#       stop("The parameters anchors and noanchors \n", "can't be positive at the same time!")
#     }
#     if (!params[["anchors"]] && !params[["noanchors"]]) {
#       stop("The parameters anchors and noanchors \n", "can't be negative at the same time!")
#     }
#   }
#   paramsCopy[["noanchors"]] <- NULL
#   params[["nocore"]] <- checkLogicalParams(
#     "nocore", params,
#     FALSE
#   )
#   if (!is.null(params[["core"]])) {
#     if (params[["core"]] && params[["nocore"]]) {
#       stop("The parameters core and nocore \n", "can't be positive at the same time!")
#     }
#     if (!params[["core"]] && !params[["nocore"]]) {
#       stop("The parameters core and nocore \n", "can't be negative at the same time!")
#     }
#   }
#   paramsCopy[["nocore"]] <- NULL
#   params[["profile"]] <- checkLogicalParams(
#     "profile", params,
#     FALSE
#   )
#   if (params[["profile"]]) {
#     if (is.null(params[["in1"]]) || is.null(params[["in2"]])) {
#       stop(
#         "The parameter profile needs the following parameters: \n",
#         "in1 and in2!"
#       )
#     }
#   }
#   paramsCopy[["profile"]] <- NULL
#   params[["refine"]] <- checkLogicalParams(
#     "refine", params,
#     FALSE
#   )
#   paramsCopy[["refine"]] <- NULL
#   params[["refinew"]] <- checkLogicalParams(
#     "refinew", params,
#     FALSE
#   )
#   paramsCopy[["refinew"]] <- NULL
#   paramsCopy[["sp"]] <- NULL
#   paramsCopy[["spn"]] <- NULL
#   params[["spscore"]] <- checkLogicalParams(
#     "spscore", params,
#     FALSE
#   )
#   paramsCopy[["spscore"]] <- NULL
#   paramsCopy[["sv"]] <- NULL
#   params[["version"]] <- checkLogicalParams(
#     "version", params,
#     FALSE
#   )
#   paramsCopy[["version"]] <- NULL
#   if (length(paramsCopy) != 0) {
#     stop(
#       "The following parameters are not known  \n", "(or have been specified",
#       "more often than once):\n    ", paste(names(paramsCopy),
#         collapse = ", ", sep = ""
#       )
#     )
#   }
#   inputSeqNames <- names(inputSeqs)
#   names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))
#   result <- .Call("RMuscle", inputSeqs, cluster, -abs(gapOpening),
#     -abs(gapExtension), maxiters, substitutionMatrix, type,
#     verbose, params,
#     PACKAGE = "msa"
#   )
#   out <- convertAlnRows(result$msa, type)
#   if (length(inputSeqNames) > 0) {
#     if (order == "aligned") {
#       perm <- match(names(out@unmasked), names(inputSeqs))
#       names(out@unmasked) <- inputSeqNames[perm]
#     }
#     else {
#       perm <- match(names(inputSeqs), names(out@unmasked))
#       out@unmasked <- out@unmasked[perm]
#       names(out@unmasked) <- inputSeqNames
#     }
#   }
#   else {
#     names(out@unmasked) <- NULL
#   }
#   standardParams <- list(
#     gapOpening = gapOpening, gapExtension = gapExtension,
#     maxiters = maxiters, verbose = verbose
#   )
#   out@params <- c(standardParams, params)
#   out@call <- deparse(sys.call())
#   out
# }
