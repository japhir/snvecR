.onLoad <- function(libname, pkgname) {
  if (getRversion() < "4.0.0") {
    backports::import("tools", "R_user_dir")
  }
}
