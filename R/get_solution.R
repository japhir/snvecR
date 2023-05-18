#' Get an Orbital Solution
#'
#' @param orbital_solution Character vector with the name of the orbital
#'   solution to use. One of `"ZB18a"` (default) from Zeebe and Lourens (2019),
#'   or `"La11"` (not yet implemented!).
# quiet and force are documented in get_ZB18a
#' @inheritParams get_ZB18a
#' @seealso [get_ZB18a()]
#' @inherit get_ZB18a references
#' @returns `get_solution()` returns a [tibble][tibble::tibble-package] with the
#'   orbital solution input and some preprocessed new columns.
#' @examples
#' \donttest{
#' get_solution()
#' }
#' @export
get_solution <- function(orbital_solution = "ZB18a", quiet = FALSE) {
  solutions <- c("ZB18a", "La11")
  if (!orbital_solution %in% solutions) {
    cli::cli_abort(c("{.var orbital_solution} must be one of: {.or {.q {solutions}}}",
                     "x" = "You've supplied {.q {orbital_solution}}"))
  }

  if (orbital_solution == "ZB18a") {
    # read in the (new?) cached result
    dat <- get_ZB18a(quiet = quiet)
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

  # if it doesn't exist, or the user wants to override the file
  if (!file.exists(rds_path) || force) {
    if (!file.exists(csv_path) || force) {
      if (!quiet) cli::cli_alert_info("The orbital solution ZB18a has not been downloaded.")

      # default to Yes downloading if not interactive (i.e. GitHub actions)
      if (!interactive()) {
        download <- TRUE
        # default to cache
        save_cache <- TRUE
      } else {# we're interactive
        # a logical, TRUE if Yes, no if otherwise
        download <- utils::menu(c("Yes", "No"), title = "Would you like to download and process it now?") == 1L
      }

      # the user would NOT like to download and process the orbital solution
      if (!download) {
        cli::cli_abort("Cannot `get_ZB18a()` without downloading the file.")
      } else {# the user would like to download and process the orbital solution
        if (interactive()) {
          save_cache <- utils::menu(c("Yes", "No"), title = "Would you like to save the results to the cache?") == 1L
        }

        if (!quiet) cli::cli_alert_info("Reading {.file {basename(raw_path)}} from website {.url {ZB18a_url}}.")

        # read the file from the website
        ZB18a <- readr::read_table(ZB18a_url,
                                   col_names = c("t", # time in days
                                                 "aa", # semimajor axis
                                                 "ee", # eccentricity
                                                 "inc", # inclination
                                                 "lph", # long periapse
                                                 "lan", # long ascending node
                                                 "arp", # argument of periapse
                                                 "mna"), # mean anomaly
                                   skip = 3,
                                   comment = "#") #|>

        if (save_cache) {
          dir.create(cachedir, recursive = TRUE, showWarnings = FALSE)
          if (!quiet) cli::cli_alert_info("The cache directory is {.file {cachedir}}.")
          # also copy the raw file to disk, even though we've read it in using read_table directly
          curl::curl_download(ZB18a_url, destfile = raw_path)
          if (!quiet) cli::cli_alert_info("Saved {.file {basename(raw_path)}} to cache.")

          # write intermediate result to csv
          readr::write_csv(ZB18a, csv_path)
          if (!quiet) cli::cli_alert_info("Saved cleaned-up {.file {basename(csv_path)}} to cache.")
        }
      }
    } else {# if we've downloaded the file but haven't prepared the solution somehow
      ZB18a <- readr::read_csv(csv_path)
    }

    # prepare the solution (i.e. calculate helper columns)
    ZB18a <- ZB18a |>
      prepare_solution()

    if (save_cache) {
      readr::write_rds(ZB18a, rds_path)
      if (!quiet) {
        cli::cli_alert("Saved solution with helper columns {.file {basename(rds_path)}} to cache.",
                       "i" = "Future runs will read from the cache, unless you specify `force = TRUE`.")
      }
    }
  } else {# if the rds file already exist, read it from the cache
    ZB18a <- readr::read_rds(rds_path)
  }

  return(ZB18a)
}
