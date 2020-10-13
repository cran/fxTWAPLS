% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fxTWAPLS.R
\name{get_pseudo}
\alias{get_pseudo}
\title{Get geographically and climatically close sites}
\usage{
get_pseudo(dist, x, cpus = 4, test_mode = FALSE, test_it = 5)
}
\arguments{
\item{dist}{distance matrix which contains the distance from other sites.}

\item{x}{the modern climate values}

\item{cpus}{number of CPUs for simultaneous iterations to execute, check
\code{parallel::detectCores()} for available CPUs on your machine.}

\item{test_mode}{boolean flag to execute the function with a limited number
of iterations, \code{test_it}, for testing purposes only.}

\item{test_it}{number of iterations to use in the test mode}
}
\value{
the geographically and climatically close sites to each test site.
}
\description{
Get the sites which are both geographically and climatically close to the
test site, which could result in pseudo-replication and inflate the
cross-validation statistics. The output will be used in
\code{\link{cv.pr.w}}
}
\examples{
# Load modern pollen data
modern_pollen <- read.csv(system.file("extdata", 
                                      "Modern_Pollen_gdd_alpha_Tmin.csv", 
                                      package = "fxTWAPLS", 
                                      mustWork = TRUE))
point <- modern_pollen[, c("Long", "Lat")]
\donttest{
test_mode <- TRUE # It should be set to FALSE before running
dist <- fxTWAPLS::get_distance(point, 
                               cpus = 2, # Remove the following line
                               test_mode = test_mode)
pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, 
                                    modern_pollen$Tmin, 
                                    cpus = 2, # Remove the following line
                                    test_mode = test_mode)
}
}
\seealso{
\code{\link{get_distance}}
}