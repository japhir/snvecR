# load the ZB18a solution from his website on load, to limit R package size
# currently commented out
## .onLoad <- function(libname, pkgname) {
##   cachedir <- tools::R_user_dir("snvecR", which = "cache")
##   ZB18apath <- paste0(cachedir, "/ZB18a.rds")
##   ZB18a_url <- "http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat"

##   if (!file.exists(ZB18apath)) {
##     cli::cli_alert_info("The orbital solution ZB18a is not downloaded yet.")

##     if (menu(c("Yes", "No"), "Would you like to download and process it now?") == 1L) {
##       cli::cli_alert_info(c("i" = "Downloading ems-plan3.dat from website.",
##                             "i" = "Will calculate some helper columns as well."))
##       dir.create(cachedir, recursive = TRUE)
##       ZB18a <- readr::read_table(ZB18a_url,
##                                  col_names = c("t", # time in days
##                                                "aa", # semimajor axis
##                                                "ee", # eccentricity
##                                                "inc", # inclination
##                                                "lph", # long periapse
##                                                "lan", # long ascending node
##                                                "arp", # argument of periapse
##                                                "mna"), # mean anomaly
##                                  skip = 3,
##                                  comment = "#") |>
##         # calculate helper columns
##         dplyr::mutate(
##           # this function is not exported so we need to use three colons
##           lphu = snvecR:::unwrap(lph),
##           lanu = snvecR:::unwrap(lan)
##         ) |>
##         dplyr::mutate(age = -t / KY2D, .after = t) |>
##         dplyr::mutate(
##           hh = ee * sin(lph / R2D),
##           kk = ee * cos(lph / R2D),
##           pp = 2 * sin(0.5 * inc / R2D) * sin(lan / R2D),
##           qq = 2 * sin(0.5 * inc / R2D) * cos(lan / R2D),
##           cc = cos(inc / R2D),
##           dd = cos(inc / R2D / 2),
##           ## /* nn <- nvec(t): normal to orbit */
##           nnx = sin(inc / R2D) * sin(lan / R2D),
##           nny = -sin(inc / R2D) * cos(lan / R2D),
##           nnz = cos(inc / R2D)
##         ) |>
##         readr::write_rds(ZB18apath)
##     }
##   } else {
##     ZB18a <- readr::read_rds(ZB18apath)
##   }
## }
