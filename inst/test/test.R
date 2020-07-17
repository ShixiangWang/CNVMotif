x = c("AEAEB", "AEAEC", "AEAEI", "AEAFA", "AEAFI", "AEBEJ", "AEBEL",
  "AEBFA", "AEBFB", "AEBFI")

show_seq_logo(x)
show_seq_logo(x, recode = TRUE)
# A:1
# B:2
# C:3
# E:1
# F:2
# I:1
# J:2
# L:4

# The 4th should be 121423

indicator <- rep(c("S", "M", "L", "E"), 6)
names(indicator) <- LETTERS[1:24]
show_seq_logo(x, recode = TRUE, indicator = indicator)
show_seq_logo(x, recode = TRUE, indicator = indicator, method = "bits")

x2 = c("AEAEB", "AEAEC", "AEAEI")
show_seq_logo(x2, recode = TRUE)

x2 = c("AEAEB", "AEAEC", "AEAEI", "AEAFA")
show_seq_logo(x2, recode = TRUE)
show_seq_logo(x2, recode = FALSE)

x2 = c("AEAEAB", "AEAEAC", "AEAEAI", "AEAFAA")
show_seq_logo(x2, recode = TRUE)
show_seq_logo(x2, recode = FALSE)
