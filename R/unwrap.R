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
  # let's stop being smart and just build it like Richard did in C
  cv <- rep(0, length(y)) # to make them the same length

  # vectorized ## dy <- diff(y) / R2D
  cx <- 0
  for (i in 2:length(y)) {
    dy <- (y[i] - y[i - 1]) / R2D # vectorized
    if (dy > pi) {
      cx <- cx - 2. * pi
    } else if (dy < -pi) {
      cx <- cx + 2. * pi
    }
    cv[i] <- cx
  }
  # ok I'll vectorize this one...
  yu <- y + cv * R2D
}
