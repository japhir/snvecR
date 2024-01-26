#' Prepare Astronomical Solution
#'
#' Calculates helper columns from an astronomical solution input.
#'
#' @export
#' @param data A data frame with the following columns:
#'
#' * `t` Time \eqn{t} (days).
#' * `ee` Eccentricity \eqn{e} (unitless).
#' * `lph` Longitude of perihelion \eqn{\varpi} (degrees).
#' * `lan` Longitude of the ascending node \eqn{\Omega} (degrees).
#' * `inc` Inclination \eqn{I} (degrees).
#' The easiest way to get this is with [get_solution()].
# inherit quiet
#' @inheritParams get_ZB
#' @returns A [tibble][tibble::tibble-package] with the new columns added.
#' @seealso [get_solution()]
#'
#' @details
#' New columns include:
#'
#' * `lphu` Unwrapped longitude of perihelion \eqn{\varpi} (degrees without
#'   jumps).
#'
#' * `lanu` Unwrapped longitude of the ascending node \eqn{\Omega} (degrees
#'   without jumps).
#'
#' * `hh` Variable: \eqn{e\sin(\varpi)}{ee * sin(lph / R2D)}.
#'
#' * `kk` Variable: \eqn{e\cos(\varpi)}{ee * cos(lph / R2D)}.
#'
#' * `pp` Variable: \eqn{2\sin(0.5I)\sin(\Omega)}{2 * sin(0.5 * inc / R2D) * n
#'   R2D)}.
#'
#' * `qq` Variable: \eqn{2\sin(0.5I)\cos(\Omega)}{2 * sin(0.5 * inc / R2D) *
#'   n  R2D)}.
#'
#' * `cc` Helper: \eqn{\cos(I)}{cos(inc / R2D)}.
#'
#' * `dd` Helper: \eqn{\cos(I)/2}{cos(inc / R2D / 2)}.
#'
#' * `nnx`, `nny`, `nnz` The \eqn{x}, \eqn{y}, and \eqn{z}-components of the
#'   Earth's orbit unit normal vector \eqn{\vec{n}}{n}, normal to Earth's
#'   instantaneous orbital plane.
# HCI = heliocentric inertial
#  \item{npx, npy, npz}{The \eqn{x}, \eqn{y}, and \eqn{z}-components of the unit
#   normal vector \eqn{\vec{n}'}{n'}, relative to ECLIPJ2000.}
# IOP = instantaneous orbit plane
prepare_solution <- function(data, quiet = FALSE) {
  mandatory_cols <- c("t", "ee", "inc", "lph", "lan")
  if (!all(mandatory_cols %in% colnames(data))) {
    if (all(c("t", "a", "e", "i", "om", "oom", "vpi", "mn") %in% colnames(data))) {
      # we're dealing with orbitN output
      cli::cli_alert_info("Renaming astronomical solution columns from orbitN syntax to snvec syntax.")
      data <- data |> dplyr::rename(ee = .data$e,
                                    inc = .data$i,
                                    lph = .data$vpi,
                                    lan = .data$oom)
    } else {
      cli::cli_abort(c(
        "Column{?s} {.col {mandatory_cols}} must be present in 'data'.",
        "x" = "Column{?s} {.col {mandatory_cols[!mandatory_cols %in% colnames(data)]}} {?is/are} missing."
      ))
    }
  }

  if (!quiet) cli::cli_alert_info("Calculating helper columns.")
  data |>
    dplyr::mutate(
      lphu = unwrap(.data$lph),
      lanu = unwrap(.data$lan)
    ) |>
    dplyr::mutate(t_ka = .data$t / KY2D, .after = "t") |>
    dplyr::mutate(age = -.data$t_ka, .after = "t_ka") |>
    dplyr::mutate(
      hh = .data$ee * sin(.data$lph / R2D),
      kk = .data$ee * cos(.data$lph / R2D),
      pp = 2 * sin(0.5 * .data$inc / R2D) * sin(.data$lan / R2D),
      qq = 2 * sin(0.5 * .data$inc / R2D) * cos(.data$lan / R2D),
      cc = cos(.data$inc / R2D),
      dd = cos(.data$inc / R2D / 2),
      ## /* nn <- nvec(t): normal to orbit */
      nnx = sin(.data$inc / R2D) * sin(.data$lan / R2D),
      nny = -sin(.data$inc / R2D) * cos(.data$lan / R2D),
      nnz = cos(.data$inc / R2D)
    )
}
