
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigminer.helper

<!-- badges: start -->

![R-CMD-check](https://github.com/ShixiangWang/sigminer.helper/workflows/R-CMD-check/badge.svg)
[![Codecov test
coverage](https://codecov.io/gh/ShixiangWang/sigminer.helper/branch/master/graph/badge.svg)](https://codecov.io/gh/ShixiangWang/sigminer.helper?branch=master)
<!-- badges: end -->

The goal of sigminer.helper is to …

## Installation

~~You can install the released version of {sigminer.helper} from
[CRAN](https://CRAN.R-project.org) with:~~

``` r
install.packages("sigminer.helper")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ShixiangWang/sigminer.helper")
```

## Example

``` r
library(sigminer.helper)
library(ggseqlogo)
#> Warning: package 'ggseqlogo' was built under R version 4.0.2
data(ggseqlogo_sample)

## Same as ggseqlogo()
p1 <- ggseqlogo2(seqs_dna[[1]])
p1
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r

## Extra feature
idor <- as.character(1:4)
names(idor) <- c("A", "C", "G", "T")
p2 <- ggseqlogo2(seqs_dna[[1]], idor = idor)
p2
```

<img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />
