# document the astronomical solutions
# we now cache them locally, so I document them with @name and return type NULL.

#' Full Astronomical Solution ZB18a for the past 100 Myr
#'
#' The HNBody output of Zeebe & Lourens (2019) after some pre-processing using
#' [prepare_solution()]. The wikipedia page on [Orbital
#' elements](https://en.wikipedia.org/wiki/Orbital_elements) describes what the
#' components relate to in order to uniquely specify an orbital plane. The
#' function asks to download the files to the user's cache directory so that
#' they can be accessed more quickly in the future.
#'
#' @format ## `get_solution("PT-ZB18a")`
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
#'   The following columns were computed from the above input with [prepare_solution()]:
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
#' * All astronomical solutions by Zeebe can be found on
#'   <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html>.
#'   They can be loaded into R quickly, using [get_solution()].
#'
#' @seealso [prepare_solution()]
#'
#' @references
#' Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
#'   Paleocene–Eocene boundary age constrained by geology and astronomy.
#'   _Science_, 365(6456), 926–929. \doi{10.1126/science.aax0612}.
#'
#' Zeebe, R. E. and Lourens, L. J. (2022). A deep-time dating tool for
#'   paleo-applications utilizing obliquity and precession cycles: The role of
#'   dynamical ellipticity and tidal dissipation. _Paleoceanography and
#'   Paleoclimatology_. \doi{10.1029/2021PA004349}
#' @name PT_ZB18a
#' @aliases PT-ZB18a
NULL

#' Astronomical Solutions ZB17 for the past 100 Myr
#'
#' @format ## `get_solution("ZB17x")`
#' A data frame with 62,501 rows and 3 columns:
#' \describe{
#'   \item{age}{Age in thousands of years before present (ka).}
#'   \item{ecc}{Eccentricity \eqn{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I} (degrees).}
#' }
#' @inherit PT_ZB18a source
#' @references
#' Zeebe, R. E. (2017). Numerical Solutions for the orbital motion of the Solar
#'   System over the Past 100 Myr: Limits and new results. _The Astronomical
#'   Journal_. \doi{10.3847/1538-3881/aa8cce}
#' @name ZB17
#' @aliases ZB17a ZB17b ZB17c ZB17d ZB17e ZB17f ZB17h ZB17i ZB17j ZB17k ZB17p
NULL

#' Astronomical Solution ZB18a for the Past 100 Myr
#'
#' @format ## `get_solution("ZB18a-100")`
#' A data frame with 62,501 rows and 3 columns:
#' \describe{
#'   \item{age}{Age in thousands of years before present (ka).}
#'   \item{ecc}{Eccentricity \eqn{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I} (degrees).}
#' }
#' @inherit PT_ZB18a source
#' @references
#' Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
#'   Paleocene–Eocene boundary age constrained by geology and astronomy.
#'   _Science_, 365(6456), 926–929. \doi{10.1126/science.aax0612}.'
#' @name ZB18a-100
NULL

#' Astronomical Solution ZB18a for the Past 300 Myr
#'
#' @format ## `get_solution("ZB18a-300")`
#' A data frame with 187,501 rows and 3 columns:
#' \describe{
#'   \item{age}{Age in thousands of years before present (ka).}
#'   \item{ecc}{Eccentricity \eqn{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I} (degrees).}
#' }
#' @inherit PT_ZB18a source
#' @references
#' Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
#'   Paleocene–Eocene boundary age constrained by geology and astronomy.
#'   _Science_, 365(6456), 926–929. \doi{10.1126/science.aax0612}.'
#' @name ZB18a-300
#' @aliases ZB18a
NULL

#' Astronomical Solutions ZB20 for the past 300 Myr
#'
#' @format ## `get_solution("ZB20x")`
#' A data frame with 187,501 rows and 3 columns:
#' \describe{
#'   \item{age}{Age in thousands of years before present (ka).}
#'   \item{ee}{Eccentricity \eqn{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I} (degrees).}
#' }
#' @inherit PT_ZB18a source
#' @references
#' Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
#'   astronomical solutions for the Cenozoic era. _Earth and Planetary Science
#'   Letters_. \doi{10.1016/j.epsl.2022.117595}
#' @name ZB20
#' @aliases ZB20a ZB20b ZB20c ZB20d
NULL
