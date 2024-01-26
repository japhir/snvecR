#' Get an Astronomical Solution
#'
# astronomical_solution, quiet, and force are documented in get_ZB
#' @inheritParams get_ZB
#' @seealso [get_ZB()]
#' @inherit ZB18a references
#' @returns A [tibble][tibble::tibble-package] with the astronomical solution
#'   (and some preprocessed new columns).
#' @examples
#' \donttest{
#' get_solution()
#' }
#' @export
get_solution <- function(astronomical_solution = "PT-ZB18a", quiet = FALSE, force = FALSE) {
  if ("data.frame" %in% class(astronomical_solution)) {
    return(prepare_solution(astronomical_solution, quiet = quiet))
  }

  # effectively this is a wrapper that checks only for valid solutions...
  solutions <- c("PT-ZB18a", # the default
                 "PT-La11", # just to thow an error
                 "ZB17a", "ZB17b", "ZB17c", "ZB17d",
                 "ZB17e", "ZB17f", "ZB17h", "ZB17i",
                 "ZB17j", "ZB17k", "ZB17p",
                 "ZB18a-100",
                 "ZB18a-300",
                 "ZB20a", "ZB20b", "ZB20c", "ZB20d")

  if (!astronomical_solution %in% solutions) {
    cli::cli_abort(c("{.var astronomical_solution} must be one of: {.or {.q {solutions}}}",
                     "x" = "You've supplied {.q {astronomical_solution}}"))
  }

  if (astronomical_solution == "PT-La11") {
    cli::cli_abort(c("i" = "Astronomical solution: La11 currently not supported.",
                     "!" = "The input OS for snvec must be either in the:",
                     "*" = "Heliocentric inertial reference frame (HCI)",
                     "*" = "Ecliptic reference frame (J2000).",
                     "x" = "The La11 solution is in the invariant/inertial reference frame.",
                     "i" = "To resolve this, you need the positions/velocities and masses of all the bodies.",
                     "i" = "Or the angles between their inertial reference frame and J2000.",
                     "i" = "Pull requests welcome."))
  }


  if (astronomical_solution == "PT-ZB18a" ||
        stringr::str_detect(astronomical_solution, "^ZB[0-9][0-9][a-z]")) {
    dat <- get_ZB(astronomical_solution, quiet = quiet, force = force)
  }

  return(dat)
}

#' Get a ZB solution
#'
#' @param astronomical_solution Character vector with the name of the desired
#'   solution.
#' @param quiet Be quiet?
#'
#'   * If `TRUE`, hide info messages.
#'
#'   * If `FALSE` (the default) print info messages and timing.
#' @param force Force re-downloading the results, even if the solution is saved
#'   to the cache.
#' @returns A [tibble][tibble::tibble-package] with the astronomical solution
#'   input and, in the case of the PT-ZB18a, some preprocessed new columns.
#' @seealso [prepare_solution()] Processes precession-tilt astronomical
#'   solution input to include helper columns.
#' @examples
#' \donttest{
#' get_ZBa()
#' }
get_ZB <- function(astronomical_solution = "PT-ZB18a",
                   quiet = FALSE,
                   force = FALSE) {
  base_url <- "http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/"
  cachedir <- tools::R_user_dir("snvecR", which = "cache")

  if (astronomical_solution == "PT-ZB18a") {
    url <- paste0(base_url, "PrecTilt/OS/ZB18a/ems-plan3.dat")
    # TODO: rename to PrecTilt, same as on ZB's website
    raw_col_names <- c("t", # time in days
                       "aa", # semimajor axis
                       "ee", # eccentricity
                       "inc", # inclination
                       "lph", # long periapse
                       "lan", # long ascending node
                       "arp", # argument of periapse
                       "mna") # mean anomaly
  } else {
    raw_col_names <- c("age", # time in ka, negative
                       "ecc", # eccentricity
                       "inc") # inclination
  }
  # different URLs for different solutions
  if (stringr::str_detect(astronomical_solution, "^ZB1[78][a-z](-100)?$")) {
    # 17 and 18 are stored in the root directory
    url <- glue::glue("{base_url}{stringr::str_replace(astronomical_solution, '-100', '')}.dat")
  } else if (stringr::str_detect(astronomical_solution, "^ZB[0-9][0-9][a-z](-300)?$")){
    # ZB20 is stored in 300Myr subdirectory
    url <- glue::glue("{base_url}300Myr/{stringr::str_replace(astronomical_solution, '-300', '')}.dat")
  }

  # where will we save our cached results?
  raw_path <- paste0(cachedir, "/", astronomical_solution, ".dat")
  csv_path <- paste0(cachedir, "/", astronomical_solution, ".csv")
  rds_path <- paste0(cachedir, "/", astronomical_solution, ".rds")

  # read final processed file from cache if available
  if (!force && file.exists(rds_path)) {
    rds <- readr::read_rds(rds_path)
    return(rds)
  }

  # read raw intermediate steps from disk
  if (!force && (file.exists(csv_path) || file.exists(raw_path))) {
    if (file.exists(csv_path)) {
      raw <- readr::read_csv(csv_path, show_col_types = FALSE)
    }

    if (file.exists(raw_path)) {
      if (astronomical_solution == "PT-ZB18a") {
        raw <- readr::read_table(raw_path,
                                 col_names = raw_col_names,
                                 skip = 3, comment = "#",
                                 show_col_types = FALSE)
      } else {
        raw <- readr::read_table(raw_path,
                                 col_names = raw_col_names,
                                 comment = "%",
                                 show_col_types = FALSE)
      }
    }
  } else {# files don't exist or force
    if (!quiet) cli::cli_alert_info("The astronomical solution {astronomical_solution} has not been downloaded.")

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

    # the user would NOT like to download and process the astronomical solution
    if (!download) {
      cli::cli_abort("Cannot `get_ZB()` without downloading the file.")
    }

    # read the file from the website
    if (!quiet) cli::cli_alert_info("Reading {.file {basename(raw_path)}} from website {.url {url}}.")
    if (astronomical_solution == "PT-ZB18a") {
      raw <- readr::read_table(url,
                               col_names = raw_col_names,
                               skip = 3, comment = "#",
                               show_col_types = FALSE)
      # calculate helper columns
      raw <- raw |> prepare_solution(quiet = quiet)
    } else {
      raw <- readr::read_table(url,
                               comment = "%",
                               col_names = raw_col_names,
                               show_col_types = FALSE)
    }

    if (!save_cache) {
      return(raw)
    } else {
      if (!dir.exists(cachedir)) {
        dir.create(cachedir, recursive = TRUE, showWarnings = TRUE)
      }
      if (!quiet) cli::cli_alert_info("The cache directory is {.file {cachedir}}.")

      # also copy the raw file to disk
      # even though we've read it in using read_table directly
      curl::curl_download(url, destfile = raw_path)
      if (!quiet) cli::cli_alert_info("Saved {.file {basename(raw_path)}} to cache.")

      # write intermediate result to csv
      readr::write_csv(raw, csv_path)
      if (!quiet) cli::cli_alert_info("Saved cleaned-up {.file {basename(csv_path)}} to cache.")

      # write final result to rds cache
      readr::write_rds(raw, rds_path)
      if (!quiet) {
        cli::cli_alert("Saved solution with helper columns {.file {basename(rds_path)}} to cache.",
                       "i" = "Future runs will read from the cache, unless you specify `force = TRUE`.")
      }
    }
    return(raw)
  }
}
