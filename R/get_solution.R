#' Get an Orbital Solution
#'
#' @param orbital_solution Character vector with the name of the orbital
#'   solution to use or a `data.frame` with a custom orbital solution.
#'   One of `"ZB18a"` (default) from Zeebe and Lourens (2019),
#'   or `"La11"` (not yet implemented!).
# quiet and force are documented in get_ZB18a
#' @inheritParams get_ZB18a
#' @seealso [get_ZB18a()]
#' @inherit get_ZB18a references
#' @returns `get_solution()` returns a [tibble][tibble::tibble-package] with
#'   the orbital solution input and some preprocessed new columns.
#' @examples
#' \donttest{
#' get_solution()
#' }
#' @export
get_solution <- function(orbital_solution = "ZB18a", quiet = FALSE, force = FALSE) {
  if ("data.frame" %in% class(orbital_solution)) {
    return(prepare_solution(orbital_solution, quiet = quiet))
  }

  solutions <- c("ZB18a", "La11")
  if (!orbital_solution %in% solutions) {
    cli::cli_abort(c("{.var orbital_solution} must be one of: {.or {.q {solutions}}}",
                     "x" = "You've supplied {.q {orbital_solution}}"))
  }

  if (orbital_solution == "ZB18a") {
    # read in the (new?) cached result
    dat <- get_ZB18a(quiet = quiet, force = force)
  }
  if (orbital_solution == "La11") {
    ## dat <- snvecR::La11
    cli::cli_abort(c("i" = "Orbital solution: La11 currently not supported.",
                     "!" = "The input OS for snvec must be either in the:",
                     "*" = "Heliocentric inertial reference frame (HCI)",
                     "*" = "Ecliptic reference frame (J2000).",
                     "x" = "The La11 solution is in the invariant/inertial reference frame.",
                     "i" = "To resolve this, you need the positions/velocities and masses of all the bodies.",
                     "i" = "Or the angles between their inertial reference frame and J2000.",
                     "i" = "Pull requests welcome."))
  }

  return(dat)
}

#' Get the ZB18a solution
#'
#' @param quiet Be quiet?
#'
#'   * If `TRUE`, hide info messages.
#'
#'   * If `FALSE` (the default) print info messages and timing.
#' @param force Force re-downloading the results, even if the solution is saved
#'   to the cache.
#' @returns `get_ZB18a()` returns a [tibble][tibble::tibble-package] with the
#'   orbital solution input and some preprocessed new columns.
#' @rdname ZB18a
#' @seealso [prepare_solution()] Processes orbital solution input to include
#'   helper columns.
#' @examples
#' \donttest{
#' get_ZB18a()
#' }
#' @export
get_ZB18a <- function(quiet = FALSE, force = FALSE) {
  ZB18a_url <- "http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat"

  cachedir <- tools::R_user_dir("snvecR", which = "cache")

  raw_path <- paste0(cachedir, "/ems-plan3.dat")
  csv_path <- paste0(cachedir, "/ZB18a.csv")
  rds_path <- paste0(cachedir, "/ZB18a.rds")

  raw_col_names <- c("t", # time in days
                     "aa", # semimajor axis
                     "ee", # eccentricity
                     "inc", # inclination
                     "lph", # long periapse
                     "lan", # long ascending node
                     "arp", # argument of periapse
                     "mna") # mean anomaly

  # read final processed file from cache if available
  if (!force && file.exists(rds_path)) {
    ZB18a <- readr::read_rds(rds_path)
    return(ZB18a)
  }

  # read raw intermediate steps from disk
  if (!force && (file.exists(csv_path) || file.exists(raw_path))) {
    if (file.exists(csv_path)) {
      ZB18a_raw <- readr::read_csv(csv_path, show_col_types = FALSE)
    }

    if (file.exists(raw_path)) {
      ZB18a_raw <- readr::read_table(raw_path,
                                     col_names = raw_col_names,
                                     skip = 3, comment = "#",
                                     show_col_types = FALSE)
    }
  } else {# files don't exist or force
    if (!quiet) cli::cli_alert_info("The orbital solution ZB18a has not been downloaded.")

    # default to downloading/caching if not interactive (i.e. GitHub actions)
    if (force || !interactive()) {
      download <- TRUE
      save_cache <- TRUE
    } else {
      # a logical, TRUE if Yes, no if otherwise
      download <- utils::menu(c("Yes", "No"),
                              title = "Would you like to download and process it now?") == 1L
      if (download) {
        save_cache <- utils::menu(c("Yes", "No"),
                                  title = "Would you like to save the results to .csv and .rds?") == 1L
      } else {
        save_cache <- FALSE
      }
    }

    # the user would NOT like to download and process the orbital solution
    if (!download) {
      cli::cli_abort("Cannot `get_ZB18a()` without downloading the file.")
    }

    # read the file from the website
    ZB18a_raw <- readr::read_table(ZB18a_url,
                                   col_names = raw_col_names,
                                   skip = 3, comment = "#",
                                   show_col_types = FALSE)
    if (!quiet) cli::cli_alert_info("Read {.file {basename(raw_path)}} from website {.url {ZB18a_url}}.")

    # calculate helper columns
    ZB18a <- ZB18a_raw |> prepare_solution(quiet = quiet)

    if (!save_cache) {
      return(ZB18a)
    } else {
      if (!dir.exists(cachedir)) {
        dir.create(cachedir, recursive = TRUE, showWarnings = TRUE)
      }
      if (!quiet) cli::cli_alert_info("The cache directory is {.file {cachedir}}.")

      # also copy the raw file to disk
      # even though we've read it in using read_table directly
      curl::curl_download(ZB18a_url, destfile = raw_path)
      if (!quiet) cli::cli_alert_info("Saved {.file {basename(raw_path)}} to cache.")

      # write intermediate result to csv
      readr::write_csv(ZB18a_raw, csv_path)
      if (!quiet) cli::cli_alert_info("Saved cleaned-up {.file {basename(csv_path)}} to cache.")

      # write final result to rds cache
      readr::write_rds(ZB18a, rds_path)
      if (!quiet) {
        cli::cli_alert("Saved solution with helper columns {.file {basename(rds_path)}} to cache.",
                       "i" = "Future runs will read from the cache, unless you specify `force = TRUE`.")
      }
    }
    return(ZB18a)
  }
}
