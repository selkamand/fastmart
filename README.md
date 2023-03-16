
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastmart

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fastmart)](https://CRAN.R-project.org/package=fastmart)
<!-- badges: end -->

The goal of fastmart is to speed up a subset of biomaRt queries by
preparing local, properly indexed sqlite databases.

## Installation

You can install the development version of fastmart like so:

``` r
install.packages("fastmart", type = "source", repos = "https://github.com/selkamand/fastmart")
```

## Quick start

The first time you use `fastmart`, two commands you need to run ONLY the
first time you use fastmart

``` r
library(fastmart)

# Initialise fastmart (creates a ~/.fastmart directory key databases will be stored)
fastmart_init()

# Cache the key tables
fastmart_cache_tables(GRCh = "37") 

# Convert HGNC to Ensemble IDs
fastmart_hgnc_to_ensembl('TBCE', chrom='1', start='235530675', end='235612283', GRCh = "37")
```
