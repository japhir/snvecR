% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{ZB20}
\alias{ZB20}
\alias{ZB20a}
\alias{ZB20b}
\alias{ZB20c}
\alias{ZB20d}
\title{Astronomical Solutions ZB20 for the past 300 Myr}
\format{
\subsection{\code{get_solution("ZB20x")}}{

A data frame with 187,501 rows and 3 columns:
\describe{
\item{time}{Time in thousands of years (kyr).}
\item{ee}{Eccentricity \eqn{e} (unitless).}
\item{inc}{Inclination \eqn{I} (degrees).}
}
}
}
\source{
All astronomical solutions by Zeebe can be found on
\url{http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}.

They can be loaded into R quickly, using \code{\link[=get_solution]{get_solution()}}.
}
\description{
Astronomical Solutions ZB20 for the past 300 Myr
}
\references{
Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
astronomical solutions for the Cenozoic era. \emph{Earth and Planetary Science
Letters}. \doi{10.1016/j.epsl.2022.117595}
}
