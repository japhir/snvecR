## unwrap function
## :PROPERTIES:
## :header-args:R: :tangle R/unwrap.R :comments org :session *R:snvec-R* :exports both :results output :eval no-export
## :END:
## :LOGBOOK:
## - State "SOME"       from              [2023-03-24 Fri 14:38]
## :END:
## [[file:snvec-3.7.5/snvec-3.7.5.c::=== unwrap()][unwrap()]]

#' unwrap angle.
#'
#' Unwrap angle. Maps jumps greater than pi to their 2pi complement.
#'
#' @param y Input vector in degrees.
#' @returns Unwrapped vector in degrees.
#' @noRd
unwrap <- function(y) {

  y <- ZB18a$lan[1:6]

  tibble::tibble(y = y) |>
    dplyr::mutate(dy = y - dplyr::lead(y, default = 0),
                  cv = 0,
           dz = ifelse(dy > 180, cv - 360,
                       ifelse(dy < -180, vc + 360, 0)),
           yu = y + cumsum(dz)) |>
    pull(yu)

  dy <- c(0, diff(y))
  dz <- dy
  dz[dy > 180] <- dy[dy > 180] - 360
  dz[dy < -180] <- dy[dy < -180] + 360
  yu <- y + cumsum(dy)
  yu
}
