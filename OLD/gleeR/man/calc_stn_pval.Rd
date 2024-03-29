% Generated by roxygen2 (4.0.2): do not edit by hand
\name{calc_stn_pval}
\alias{calc_stn_pval}
\title{calc_stn_pval.}
\usage{
calc_stn_pval(A, B, m, num_iter)
}
\arguments{
\item{A}{Replicates in the first condition.}

\item{B}{Replicates in the second condition.}

\item{m}{List containing model-fit information, as generated by \code{\link{fit_model}}}

\item{num_iter}{Number of iterations for generating STN distribution}
}
\value{
List containing model-based STN and p-value for each protein
}
\description{
Calculate model-based signal-to-noise (STN) and p-values.
}

