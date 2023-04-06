## document the data

#' Astronomical Solution from Zeebe et al., 2018
#'
#' The HNBody output of Zeebe & Lourens 2019 \url{10.1126/science.aax0612}.
#'
#' The wikipedia page on [Orbital elements](https://en.wikipedia.org/wiki/Orbital_elements)
#' describes what the components relate to to uniquely specify an orbital plane.
#'
#' @format ## `ZB18a`
#' A data frame with 250,001 rows and 8 columns:
#' \describe{
#'   \item{t}{Time in days.}
#'   \item{age}{Age in thousands of years (kyr) before present.}
#'   \item{aa}{Semimajor axis.}
#'   \item{ee}{Eccentricity}
#'   \item{inc}{Inclination}
#'   \item{lph}{Long periapse}
#'   \item{lan}{Long ascending node}
#'   \item{arp}{Argument of periapse}
#'   \item{mna}{Mean anomaly}
#'   \item{hh}{helper hh}
#'   \item{kk}{helper hk}
#'   \item{pp}{helper pp}
#'   \item{qq}{helper qq}
#'   \item{cc}{helper cc}
#'   \item{dd}{helper dd}
#'   \item{nnx, nny, nnz}{Vector of Earth's orbit normal.}
#' }
#' @source <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html>
"ZB18a"
