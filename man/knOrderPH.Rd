% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/k_of_n_order_ph.R
\name{knOrderPH}
\alias{knOrderPH}
\title{Generate the sub-transition matrix of the k:n order statistic given \eqn{X_1,..., X_n} phase-types.}
\usage{
knOrderPH(listAlpha, listS, k)
}
\arguments{
\item{listAlpha}{list of initial distributions of the n phase-type distributions.}

\item{listS}{list of sub-intensity matrices of the n phase-type distributions.}

\item{k}{the order statistic of \eqn{n}.}
}
\value{
Sub-intensity matrix and initial vector of phase-type representation of kth order statistic. Also dimensions of each block matrix.
}
\description{
Generate the sub-transition matrix of the k:n order statistic given \eqn{X_1,..., X_n} phase-types.
}
