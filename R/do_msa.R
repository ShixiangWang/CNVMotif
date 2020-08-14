#' Run Modified Multiple Sequence Alignment
#'
#' @param x a character vector.
#' @param substitutionMatrix substitution matrix for scoring matches and mismatches.
#' Default is `NULL`, use `1` for match and `0` for unmatch.
#' @param gapOpening gap opening penalty; Note that the sign of this parameter is ignored.
#' @param gapExtension gap extension penalty; Note that the sign of this parameter is ignored.
#' @param verbose if `TRUE`, print extra info.
#' @param ... other arguments passing to [msa::msa]
#'
#' @return a `list`.
#' @export
#'
#' @examples
#' r <- do_msa(c("ABCDF", "BCDEF"))
#' r
#' @testexamples
#' expect_is(r, "list")
do_msa <- function(x, substitutionMatrix = NULL,
                   gapOpening = 6,
                   gapExtension = 1,
                   verbose = FALSE, ...) {

  if (is.null(substitutionMatrix)) {
    substitutionMatrix <- matrix(0, nrow = 6, ncol = 6,
                                 dimnames = list(LETTERS[1:6], LETTERS[1:6]))
    diag(substitutionMatrix) <- 1
  }
  s <- Biostrings::AAStringSet(x)
  names(s) <- paste0("s", seq_along(x))
  msa <- msaMuscleCustom(s,
                         gapOpening = gapOpening,
                         gapExtension = gapExtension,
                         substitutionMatrix = substitutionMatrix,
                         verbose = verbose,
                         type = "custom",
                         ...)

  ConsensusSequence <- msa::msaConsensusSequence(msa)
  ConservationScore <- msa::msaConservationScore(msa, substitutionMatrix)
  Seqs <- Biostrings::unmasked(msa)
  list(
    MSA = Seqs,
    ConsensusSequence = ConsensusSequence,
    ConservationScore = ConservationScore
  )
}
