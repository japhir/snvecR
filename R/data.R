## document the data

#' Orbital Solution from Zeebe et al., 2018.
#'
#' The HNBody output of Zeebe & Lourens 2019.
#'
#' The wikipedia page on [Orbital elements](https://en.wikipedia.org/wiki/Orbital_elements)
#' describes what the components relate to in order to uniquely specify an orbital plane.
#'
#' @format ## `ZB18a`
#' A data frame with 250,001 rows and 8 columns:
#' \describe{
#'   \item{t}{Time in days.}
#'   \item{age}{Age in thousands of years (kyr) before present.}
#'   \item{aa}{Semimajor axis.}
#'   \item{ee}{Eccentricity.}
#'   \item{inc}{Inclination.}
#'   \item{lph}{Long periapse.}
#'   \item{lan}{Long ascending node.}
#'   \item{arp}{Argument of periapse.}
#'   \item{mna}{Mean anomaly.}
#'   \item{hh}{helper hh.}
#'   \item{kk}{helper hk.}
#'   \item{pp}{helper pp.}
#'   \item{qq}{helper qq.}
#'   \item{cc}{helper cc.}
#'   \item{dd}{helper dd.}
#'   \item{nnx, nny, nnz}{Vector of Earth's orbit normal.}
#' }
#' @source
#' All orbital solutions can be found on <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html>.
#' The specific one we use here is available at <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat>.
#'
#' @references
#' Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
#'   Paleocene–Eocene boundary age constrained by geology and astronomy.
#'   _Science_, 365(6456), 926–929. <https://doi.org/10.1126/science.aax0612>
"ZB18a"
