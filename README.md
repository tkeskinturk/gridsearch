
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gridsearch <img src="sticker/gridsearch.png" align="right" alt="" width="120" />

<!-- badges: start -->
<!-- badges: end -->

`gridsearch` is an R package used for adjudicating plausible data
generation processes (DGP) for a given panel dataset. In its current
version, `gridsearch` allows researchers to:

- Generate simulated datasets from a given set of DGP parameters,
- Implement a simple grid-search algorithm to see what plausible
  combinations of parameter values may approximate to the observed
  marginal distributions,
- Plot distance values generated from the grid-search algorithm,
- Generate simple simulations for calculating classification accuracy
  measures to see whether, for a given DGP, we can identify individual
  changers with a certain level of accuracy.

This package is a companion to the article, *The Promises and Pitfalls
of Using Panel Data to Understand Individual Belief Change*, stored in
[SocArXiv](https://osf.io/preprints/socarxiv/rhf4q).

## Installation

You can install the development version of `gridsearch` from Github
with:

``` r
# install.packages("devtools")
devtools::install_github("tkeskinturk/gridsearch")
```
