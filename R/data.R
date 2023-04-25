## document the data

#' Orbital Solution ZB18a
#'
#' The HNBody output of Zeebe & Lourens (2019).
#'
#' The wikipedia page on [Orbital elements](https://en.wikipedia.org/wiki/Orbital_elements)
#' describes what the components relate to in order to uniquely specify an orbital plane.
#'
#' @format ## `ZB18a`
#' A data frame with 250,001 rows and 20 columns:
#' \describe{
#'   \item{t}{Time \eqn{t} (days).}
#'   \item{age}{Age in thousands of years before present (ka).}
#'   \item{aa}{Semimajor axis \eqn{a} in astronomical units (au).}
#'   \item{ee}{Eccentricity \eqn{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I} (degrees).}
#'   \item{lph}{Longitude of perihelion \eqn{\varpi} (degrees).}
#'   \item{lan}{Longitude of the ascending node \eqn{\Omega} (degrees).}
#'   \item{arp}{Argument of perihelion \eqn{\omega} (degrees).}
#'   \item{mna}{Mean anomaly \eqn{M} (degrees).}
#'
#'   The following columns were computed from the above input:
#'
#'   \item{lphu}{Unwrapped longitude of perihelion \eqn{\varpi} (degrees without
#'   jumps).}
#'
#'   \item{lanu}{Unwrapped longitude of the ascending node \eqn{\Omega}
#'   (degrees without jumps).}
#'
#'   \item{hh}{Variable: \eqn{e\sin(\varpi)}{ee * sin(lph / R2D)}.}
#'
#'   \item{kk}{Variable: \eqn{e\cos(\varpi)}{ee * cos(lph / R2D)}.}
#'
#'   \item{pp}{Variable: \eqn{2\sin(0.5I)\sin(\Omega)}{2 * sin(0.5 * inc / R2D) *
#'   sin(lan / R2D)}.}
#'
#'   \item{qq}{Variable: \eqn{2\sin(0.5I)\cos(\Omega)}{2 * sin(0.5 * inc / R2D) *
#'   cos(lan / R2D)}.}
#'
#'   \item{cc}{Helper: \eqn{\cos(I)}{cos(inc / R2D)}.}
#'
#'   \item{dd}{Helper: \eqn{\cos(I)/2}{cos(inc / R2D / 2)}.}
#'
#'   \item{nnx, nny, nnz}{The \eqn{x}, \eqn{y}, and \eqn{z}-components of the
#'   Eart's orbit unit normal vector \eqn{\vec{n}}{n}, normal to Earth's
#'   instantaneous orbital plane.}
#    HCI = heliocentric inertial
#   \item{npx, npy, npz}{The \eqn{x}, \eqn{y}, and \eqn{z}-components of the
#   unit normal vector \eqn{\vec{n}'}{n'}, relative to ECLIPJ2000.}
#   IOP = instantaneous orbit plane
#'
#' }
#' @source
#'
#' * All orbital solutions by Zeebe can be found on
#'   <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html>.
#'
#' * The specific one we use here is available at
#' <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat>.
#'
#' @references
#'
#' Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
#'   Paleocene–Eocene boundary age constrained by geology and astronomy.
#'   _Science_, 365(6456), 926–929. \doi{10.1126/science.aax0612}.
"ZB18a"
