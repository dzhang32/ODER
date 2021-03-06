% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimise_ERs.R
\name{calc_ERs_delta}
\alias{calc_ERs_delta}
\title{Calculates delta for sets of ERs}
\usage{
calc_ERs_delta(ERs_MCCs_MRGs, opt_gr, delta_func = .delta)
}
\arguments{
\item{ERs_MCCs_MRGs}{Sets of ERs across various MCCs/MRGs - output of
\code{\link{gen_ERs}}.}

\item{opt_gr}{GRanges object that contains the regions that ideally, you want
to the ER definitions to match}

\item{delta_func}{Function that calculates the delta between ERs and
\code{opt_gr}. Takes as input a set of ERs from \code{ERs_MCCs_MRGs} and
\code{opt_gr}. Then outputs a tibble/dataframe containing the summarised
delta scores for that set of one set of ERs.}
}
\value{
tibble/dataframe containing summarised delta values. One row per set
  of ERs.
}
\description{
\code{calc_ERs_delta} calculates the delta/difference between a set of ERs
and another given set of GRanges that are the optimal
}
