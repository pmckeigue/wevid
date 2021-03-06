% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wdists.R
\name{weightsofevidence}
\alias{weightsofevidence}
\title{Calculate weights of evidence in natural log units}
\usage{
weightsofevidence(posterior.p, prior.p)
}
\arguments{
\item{posterior.p}{Vector of posterior probabilities generated by using model
to predict on test data}

\item{prior.p}{Prior probabilities on test data.}
}
\value{
The weight of evidence in nats for each observation.
}
\description{
Calculate weights of evidence in natural log units
}
\examples{
data("cleveland") # load example dataset
W <- with(cleveland, weightsofevidence(posterior.p, prior.p))
densities.unadj <- Wdensities.unadjusted(cleveland$y, W)
densities.adj <- Wdensities.fromraw(densities.unadj)
plotWdists(densities.unadj, densities.adj)

}
