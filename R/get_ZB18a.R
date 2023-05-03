#' Get an Orbital Solution
#'
#' @param orbital_solution Character vector with the name of the orbital
#'   solution to use. One of `"ZB18a"` (default) from Zeebe and Lourens (2019),
#'   or `"La11"` (not yet implemented!).
get_solution <- function(orbital_solution = "ZB18a") {
  solutions <- c("ZB18a", "La11")
  if (!orbital_solution %in% solutions) {
    cli::cli_abort(c("{.var orbital_solution} must be one of: {.or {.q {solutions}}}",
                     "x" = "You've supplied {.q {orbital_solution}}"))
  }

  if (orbital_solution == "ZB18a") {
    # read in the (new?) cached result
    dat <- get_ZB18a()
  }
  if (orbital_solution == "La11") {
    ## dat <- snvecR::La11
    cli::cli_abort(c("Orbital solution: La11 currently not supported.",
      "!" = "The input OS for snvec must be in the Heliocentric Inertial Reference frame (HCI) (J2000).",
      "x" = "The La11 solution is in the invariant reference frame.",
      "i" = "Pull requests welcome."))
  }

  return(dat)
}

#' Get the ZB18a solution from the cache
#'
# The cache is created in when the package is loaded, in zzz.R
get_ZB18a <- function() {
  cachedir <- tools::R_user_dir("snvecR", which = "cache")
  ZB18apath <- paste0(cachedir, "/ZB18a.rds")

  # if you picked no during the onLoad, I assume read_rds will tell you the
  # file doesn't exist
  ZB18a <- readr::read_rds(ZB18apath)
  return(ZB18a)
}
