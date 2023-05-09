# using approxfun
# http://desolve.r-forge.r-project.org/ has an article on time-varying inputs
# we can use approxfun to generate a function that approximates =col= for timestep t.

#' Linearly interpolate a dataframe.
#'
#' @param dat The dataframe.
#' @param col A character vector with the column we want to interpolate on.
#' @returns A function that predicts `col` as a function of `t` in `dat`.
#'
#' @seealso [qinterp()] for quick interpolation of a single timestep.
#' @noRd
approxdat <- function(dat, col) {
  # I'm not putting any input checks because it's an internal function
  dat |>
    dplyr::select(all_of(c("t", col))) |>
    approxfun(rule = 2)
}

# implement qinterp similar to the C-routine
#' Quickly interpolate a single value.
#'
#' `qinterp()` is a custom linear interpolation algorithm that is much faster
#' than using the full vectorized `[approx()]` or `[approxfun()]`, because it
#' only interpolates the single value of the current timestep.
#'
#' @param y The vector to interpolate.
#' @param ds The difference in timestep in the astronomical solution.
#' @param dx The difference between the current timestep and the timestep in the astronomical solution.
#' @param m The index variable of the current position in the astronomical solution.
#' @returns The vector of interpolated results.
#' @examples
#' # interpolate ZB18a$lph[[1:4]]
#' qinterp(ZB18a$lph[[1:4]], ds = -146100, dx =  -18262.5, m = 2)
#'
#' @seealso [approxdat()] for linear interpolation of the full astronomical solution.
#' @noRd
qinterp <- function(y, ds, dx, m) {
  yi <- y[m]
  dy <- 0
  dsa <- abs(ds)
  dxa <- abs(dx)
  mm <- 1L
  ## if (dxa > 0) {
  mm <- m - as.integer(sign(dx))
  dy <- y[mm] - y[m]
  yi <- yi + dy * dxa / dsa
  ## }
  return(yi)
}
