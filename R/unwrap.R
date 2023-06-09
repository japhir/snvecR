#' Unwrap angle.
#'
#' Maps jumps greater than pi to their 2pi complement.
#'
#' @param y Input vector in degrees.
#' @returns Unwrapped vector in degrees.
#' @noRd
unwrap <- function(y) {
  ## y <- ZB18a$lan[1:6]
  cx <- 0.0 # single number, c in C but reserved namespace
  dy <- 0.0 # single number
  # vector
  cv <- rep(0, length(y))

  for (i in 2:length(y)) {
    dy <- (y[i] - y[i - 1]) / R2D
    if (dy > pi) {
      cx <- cx - 2 * pi
    } else if (dy < -pi) {
      cx <- cx + 2 * pi
    }
    cv[i] <- cx
  }
  yu <- rep(0, length(y))
  for (i in seq_along(y)) {
    yu[i] <- y[i] + cv[i] * R2D
  }
  yu
}
