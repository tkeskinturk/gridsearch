% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gridSearch.R
\name{gridSearch}
\alias{gridSearch}
\title{Grid Search for Adjudicating DGPs}
\usage{
gridSearch(
  data,
  yname,
  tname,
  pname,
  n_samples = 1000,
  ic_min = 0,
  ic_max = 2,
  pc_min = 0,
  pc_max = 1,
  bal_min = 0,
  bal_max = 1,
  rel_min = 0,
  rel_max = 1,
  fix = "none",
  fix_at = 1,
  verbose = TRUE
)
}
\arguments{
\item{data}{A panel dataframe in the long format.}

\item{yname}{The outcome identifier.}

\item{tname}{The time identifier.}

\item{pname}{The unit identifier.}

\item{n_samples}{The number of samples from priors.}

\item{ic_min}{Minimum strength of change parameter.}

\item{ic_max}{Maximum strength of change parameter.}

\item{pc_min}{Minimum rate of change parameter.}

\item{pc_max}{Maximum rate of change parameter.}

\item{bal_min}{Minimum balance parameter.}

\item{bal_max}{Maximum balance parameter.}

\item{rel_min}{Minimum reliability parameter.}

\item{rel_max}{Maximum reliability parameter.}

\item{fix}{Fix any parameter? (values: "ic_sample," "pc_sample," etc.).}

\item{fix_at}{At what value to fix it on?}

\item{verbose}{Whether to see detailed messages.}
}
\value{
A data frame.
}
\description{
This function uses an Approximate Bayesian Computation (ABC) algorithm to calculate the extent to which a dataframe can be approximated by a known DGP. It implements simulations across plausible values---rate of change, strength of change, directionality, and reliability---and provides a dataframe with an error score per DGP.
}
