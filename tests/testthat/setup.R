.onLoad <- function(libname, pkgname) {
  if (getRversion() < "4.0.0") {
    backports::import("tools", "R_user_dir")
  }
}

cleanup <- function() {
  cachedir <- R_user_dir("snvecR", which = "cache")
  if (dir.exists(cachedir) && !interactive()) {
    cli::cli_inform("Removing {.file {cachedir}} from reproducible environment.")
    unlink(cachedir, recursive = TRUE)
  }
}

withr::defer(cleanup(), teardown_env())
