utils::globalVariables(
  c("sample", "start", "end", "lenVal", "segVal", "ID", "chromosome", "%dopar%", "foreach", "Seqs")
)

# Internal functions from other packages ----------------------------------

bits_method <- getFromNamespace("bits_method", "ggseqlogo")
probability_method <- getFromNamespace("probability_method", "ggseqlogo")
get_font <- getFromNamespace("get_font", "ggseqlogo")
get_col_scheme <- getFromNamespace("get_col_scheme", "ggseqlogo")
guessSeqType <- getFromNamespace("guessSeqType", "ggseqlogo")
matrix_to_heights <- getFromNamespace("matrix_to_heights", "ggseqlogo")
newRange <- getFromNamespace("newRange", "ggseqlogo")
