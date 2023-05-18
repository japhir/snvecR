# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(snvecR)

test_check("snvecR")

# clean up after tests
cachedir <- tools::R_user_dir("snvecR", which = "cache")
if (dir.exists(cachedir)) {
  cli::cli_inform("Removing {.file {cachedir}} from reproducible environment.")
  unlink(cachedir, recursive = TRUE)
}
