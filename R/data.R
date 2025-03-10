# document the astronomical solutions
# we cache them locally, so I document them with @name and return type NULL.

#' Full Astronomical Solution ZB18a
#'
#' The full ZB18a solution spans the past 100 Myr.
#' It contains the HNBody output of Zeebe & Lourens (2019) after some pre-processing using
#' [prepare_solution()]. The wikipedia page on [Orbital
#' elements](https://en.wikipedia.org/wiki/Orbital_elements) describes what the
#' components relate to in order to uniquely specify an orbital plane.
#'
#' @format ## `get_solution("full-ZB18a")`
#' A data frame with 250,001 rows and 20 columns:
#' \describe{
#'   \item{t}{Time \eqn{t}{t} (days).}
#'   \item{time}{Time in thousands of years (kyr).}
#'   \item{aa}{Semimajor axis \eqn{a}{a} in astronomical units (au).}
#'   \item{ee}{Eccentricity \eqn{e}{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I}{I} (degrees).}
#'   \item{lph}{Longitude of perihelion \eqn{\varpi}{varpi} (degrees).}
#'   \item{lan}{Longitude of the ascending node \eqn{\Omega}{Omega} (degrees).}
#'   \item{arp}{Argument of perihelion \eqn{\omega}{omega} (degrees).}
#'   \item{mna}{Mean anomaly \eqn{M}{M} (degrees).}
#'
#'   The following columns were computed from the above input with [prepare_solution()]:
#'
#'   \item{lphu}{Unwrapped longitude of perihelion \eqn{\varpi}{varpi} (degrees
#'   without jumps).}
#'
#'   \item{lanu}{Unwrapped longitude of the ascending node \eqn{\Omega}{Omega}
#'   (degrees without jumps).}
#'
#'   \item{hh}{Variable: \eqn{e\sin(\varpi)}{e sin(lph)}.}
#'
#'   \item{kk}{Variable: \eqn{e\cos(\varpi)}{e cos(lph)}.}
#'
#'   \item{pp}{Variable: \eqn{2\sin(0.5I)\sin(\Omega)}{2 sin(0.5 I)
#'   sin(Omega)}.}
#'
#'   \item{qq}{Variable: \eqn{2\sin(0.5I)\cos(\Omega)}{2 sin(0.5 inc) *
#'   cos(Omega)}.}
#'
#'   \item{cc}{Helper: \eqn{\cos(I)}{cos(I)}.}
#'
#'   \item{dd}{Helper: \eqn{\cos(I)/2}{cos(I) / 2}.}
#'
#'   \item{nnx, nny, nnz}{The \eqn{x}{x}, \eqn{y}{y}, and \eqn{z}{z}-components of the
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
#' All astronomical solutions by Zeebe can be found on
#'   <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html>.
#'
#' They can be loaded into R quickly, using [get_solution()].
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
#' @name full_ZB18a
#' @aliases full-ZB18a
NULL

#' Astronomical Solutions ZB17
#'
#' The ZB17x eccentricity solutions span the past 100 Myr.
#' Available solutions include `"ZB17a"`, `"ZB17b"`, `"ZB17c"`, `"ZB17d"`,
#' `"ZB17e"`, `"ZB17f"`, `"ZB17h"`, `"ZB17i"`, `"ZB17j"`, `"ZB17k"`, and
#' `"ZB17p"`.
#'
#' @format ## `get_solution("ZB17x")`
#' A data frame with 62,501 rows and 3 columns:
#' \describe{
#'   \item{time}{Time in thousands of years (kyr).}
#'   \item{ecc}{Eccentricity \eqn{e}{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I}{I} (degrees).}
#' }
#' @inherit full_ZB18a source
#' @references
#' Zeebe, R. E. (2017). Numerical Solutions for the orbital motion of the Solar
#'   System over the Past 100 Myr: Limits and new results. _The Astronomical
#'   Journal_. \doi{10.3847/1538-3881/aa8cce}
#' @name ZB17
#' @aliases ZB17a ZB17b ZB17c ZB17d ZB17e ZB17f ZB17h ZB17i ZB17j ZB17k ZB17p
NULL

#' Astronomical Solution ZB18a
#'
#' The ZB18a_100 eccentricity solution spans the past 100 Myr. See [ZB18a_300]
#' for the past 300 Myr.
#'
#' @format ## `get_solution("ZB18a-100")`
#' A data frame with 62,501 rows and 3 columns:
#' \describe{
#'   \item{time}{Time in thousands of years (kyr).}
#'   \item{ecc}{Eccentricity \eqn{e}{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I}{I} (degrees).}
#' }
#' @inherit full_ZB18a source
#' @references
#' Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
#'   Paleocene–Eocene boundary age constrained by geology and astronomy.
#'   _Science_, 365(6456), 926–929. \doi{10.1126/science.aax0612}.'
#' @name ZB18a_100
#' @aliases "ZB18a-100"
NULL

#' Astronomical Solution ZB18a
#'
#' The ZB18a_300 eccentricity solution spans the past 300 Myr. See [ZB18a_100]
#' for the past 100 Myr.
#'
#' @format ## `get_solution("ZB18a-300")`
#' A data frame with 187,501 rows and 3 columns:
#' \describe{
#'     \item{time}{Time in thousands of years (kyr).}
#'     \item{ecc}{Eccentricity \eqn{e}{e} (unitless).} \item{inc}{Inclination
#'     \eqn{I}{I} (degrees).}
#' }
#' @inherit full_ZB18a source
#' @references Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and
#'   the Paleocene–Eocene boundary age constrained by geology and astronomy.
#'   _Science_, 365(6456), 926–929. \doi{10.1126/science.aax0612}.'
#'
#' Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
#'   astronomical solutions for the Cenozoic era. _Earth and Planetary Science
#'   Letters_. \doi{10.1016/j.epsl.2022.117595}
#' @name ZB18a_300
#' @aliases ZB18a "ZB18a-300"
NULL

#' Astronomical Solutions ZB20
#'
#' The ZB20x eccentricity solutions span the past 300 Myr. Available solutions
#' include `"ZB20a"`, `"ZB20b"`, `"ZB20c"`, and `"ZB20d"`.
#'
#' @format ## `get_solution("ZB20x")`
#' A data frame with 187,501 rows and 3 columns:
#' \describe{
#'   \item{time}{Time in thousands of years (kyr).}
#'   \item{ee}{Eccentricity \eqn{e}{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I}{I} (degrees).}
#' }
#' @inherit full_ZB18a source
#' @references
#' Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
#'   astronomical solutions for the Cenozoic era. _Earth and Planetary Science
#'   Letters_. \doi{10.1016/j.epsl.2022.117595}
#' @name ZB20
#' @aliases ZB20a ZB20b ZB20c ZB20d
NULL

#' Astronomical Solutions PT-ZB18a(E<sub>d</sub>,T<sub>d</sub>)
#'
#' The pre-computed precession-tilt solutions
#' PT-ZB18a(E<sub>d</sub>,T<sub>d</sub>) span the past 100 Myr. Available
#' solutions include all combinations of dynamical ellipticity values between
#' 0.9950 and 1.0050 in increments of 0.0010 and tidal dissipation values
#' between 0.000 and 1.2000 in increments of 0.1.
#'
#' @format ## `get_solution("PT-ZB18a(1.000,1.000)")`
#' A data frame with 249,480 rows and 4 columns:
#' \describe{
#'   \item{time}{Time in thousands of years (kyr).}
#'   \item{epl}{Obliqity \eqn{\epsilon}{epsilon} (radians).}
#'   \item{phi}{Axial Precession \eqn{\phi}{phi} (radians).}
#'   \item{cp}{Climatic Precession \eqn{e \sin(\bar{\omega})}{e sin(baromega)} (unitless).}
#' }
#' @inherit full_ZB18a source
#' @references
#' Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
#'   astronomical solutions for the Cenozoic era. _Earth and Planetary Science
#'   Letters_. \doi{10.1016/j.epsl.2022.117595}
#' @name PT_ZB18a
#' @aliases PT-ZB18a "PT-ZB18a(1,1)" "PT-ZB18a(1,0)"
NULL

#' Astronomical Solutions ZB23.Rxx
#'
#' The ZB23.Rxx eccentricity solutions spand the past 3.6 Gyr. Available solutions include
#' `"ZB23.R01"` to `"ZB23.R60"` and `"ZB23.R62"` to `"ZB23.R64"`.
#'
#' @format ## `get_solution("ZB23.Rxx")`
#' A data frame with 8,750,001 rows and 5 columns:
#' \describe{
#'   \item{time}{Time in thousands of years (kyr).}
#'   \item{ecc}{Eccentricity \eqn{e}{e} (unitless).}
#'   \item{inc}{Inclination \eqn{I}{I} (radians).}
#'   \item{obliquity}{Obliqity \eqn{\epsilon}{epsilon} (radians).}
#'   \item{cp}{Climatic Precession \eqn{e \sin(\bar{\omega})}{e sin(baromega)} (unitless).}
#' }
#' @inherit full_ZB18a source
#' @references
#' Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
#'   astronomical solutions for the Cenozoic era. _Earth and Planetary Science
#'   Letters_. \doi{10.1016/j.epsl.2022.117595}
#' @name ZB23
#' @aliases ZB23.Rxx
#' @aliases ZB23.R01 ZB23.R02 ZB23.R03 ZB23.R04 ZB23.R05 ZB23.R06 ZB23.R07 ZB23.R08 ZB23.R09 ZB23.R10
#' @aliases ZB23.R11 ZB23.R12 ZB23.R13 ZB23.R14 ZB23.R15 ZB23.R16 ZB23.R17 ZB23.R18 ZB23.R19 ZB23.R20
#' @aliases ZB23.R21 ZB23.R22 ZB23.R23 ZB23.R24 ZB23.R25 ZB23.R26 ZB23.R27 ZB23.R28 ZB23.R29 ZB23.R30
#' @aliases ZB23.R31 ZB23.R32 ZB23.R33 ZB23.R34 ZB23.R35 ZB23.R36 ZB23.R37 ZB23.R38 ZB23.R39 ZB23.R40
#' @aliases ZB23.R41 ZB23.R42 ZB23.R43 ZB23.R44 ZB23.R45 ZB23.R46 ZB23.R47 ZB23.R48 ZB23.R49 ZB23.R50
#' @aliases ZB23.R51 ZB23.R52 ZB23.R53 ZB23.R54 ZB23.R55 ZB23.R56 ZB23.R57 ZB23.R58 ZB23.R59 ZB23.R60
#' @aliases ZB23.R62 ZB23.R63 ZB23.R64
# note that 61 is missing!
NULL
