.onLoad <- function(libname, pkgname) {
    cachedir <- tools::R_user_dir("snvecR", which = "cache")
    csv_path <- paste0(cachedir, "/ZB18a.csv")
    rds_path <- paste0(cachedir, "/ZB18a.rds")

    # if it doesn't exist
    if (!file.exists(rds_path)) {
      if (!file.exists(csv_path)) {
        cli::cli_alert_info("The orbital solution ZB18a has not been downloaded.")
        ZB18a_url <- "http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat"

        # ask user to download it
        if (interactive()) {
          choice <- utils::menu(c("Yes", "No"), title = "Would you like to download and process it now?")
        } else {
          # default to Yes if not interactive (i.e. GitHub actions)
          choice <- 1L
        }

        if (choice == 1L) {
          dir.create(cachedir, recursive = TRUE, showWarnings = FALSE)
          cli::cli_alert_info("Downloading ems-plan3.dat from website {ZB18a_url}.")

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
                                     comment = "#") |>
            # write intermediate result to csv
            readr::write_csv(csv_path)
          cli::cli_alert_info("Saved cleaned-up file as ZB18a.csv")
        }
      } else {# if we've downloaded the file but haven't prepared the solution somehow
          ZB18a <- readr::read_csv(csv_path)
      }

      # prepare the solution and save to rds
      ZB18a <- ZB18a |>
        prepare_solution() |>
        readr::write_rds(rds_path)
      cli::cli_alert_info("Saved solution with helper columns to {rds_path}")
    }
}
