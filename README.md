
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastmart

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fastmart)](https://CRAN.R-project.org/package=fastmart)
[![R-CMD-check](https://github.com/selkamand/fastmart/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/fastmart/actions/workflows/R-CMD-check.yaml)
![r-universe](https://selkamand.r-universe.dev/badges/fastmart)

<!-- badges: end -->

The goal of fastmart is to speed up a subset of biomart queries by
caching local, properly indexed sqlite databases.

## Installation

You can install fastmart from R-Universe:

``` r
install.packages("fastmart", repos = "https://ropensci.r-universe.dev")
```

Or install the develeopment version from github:

``` r
# install.packages('remotes')
remotes::install_github("selkamand/fastmart")
```

## Quick start

The first time you use this package, you will need to **initialise**
fastmart and cache the tables we need. This can be done as follows

``` r
library(fastmart)

# Initialise fastmart (creates a ~/.fastmart directory key databases will be stored)
fastmart_init()

# Cache the key tables
fastmart_cache_tables(GRCh = "37") 
```

You can now use any of the supported quick lookups

## Available Functionality

### Convert Gene IDs

Convert HGNC Symbols to ENSEMBL gene IDs

``` r
# Convert HGNC to Ensemble IDs
fastmart_convert_hgnc_to_ensembl('TBCE', chrom='1', start='235530675', end='235612283', GRCh = "37")
  # > ENSG00000116957.8
```

### Annotate Gene with biotype

``` r
fastmart_annotate_biotype('ENSG00000282455.1')
#> IG_D_gene
```
