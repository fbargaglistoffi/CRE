% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_rules.R
\name{select_rules}
\alias{select_rules}
\title{Select rules}
\usage{
select_rules(
  rules_matrix,
  rules_list,
  ite,
  stability_selection,
  cutoff,
  pfer,
  penalty_rl
)
}
\arguments{
\item{rules_matrix}{The rules matrix.}

\item{rules_list}{A list of rules.}

\item{ite}{The estimated ITE.}

\item{stability_selection}{Whether or not using stability selection.}

\item{cutoff}{Threshold (percentage) defining the minimum cutoff value for
the stability scores. Only for stability selection.}

\item{pfer}{Upper bound for the per-family error rate (tolerated amount of
falsely selected rules). Only for stability selection.}

\item{penalty_rl}{Order of penalty for rules length during LASSO
regularization (i.e. 0: no penalty, 1: rules_length, 2: rules_length^2).}
}
\value{
A minimal set of rules linearly decomposing the CATE.
}
\description{
Given a set of rules, selects the minimal set linearly decomposing the
Conditional Average Treatment Effect (CATE) by LASSO (optionally with
Stability Selection).
}
\keyword{internal}
