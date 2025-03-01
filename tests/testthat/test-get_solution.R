test_that("get_solution() bad inputs throw errors", {
  pth <- withr::local_tempdir(pattern = "snvecR")
  withr::local_options(list(snvecR.cachedir = pth))
  expect_error(get_solution(astronomical_solution = "hoi"))
  expect_error(get_solution(astronomical_solution = "full-La11"))
  expect_error(get_solution(astronomical_solution = data.frame(x = 1, y = 5))) # this runs via prepare_solution
})

test_that("get_solution() can return a dataframe", {
  ZB18a_head <- structure(list(t = c(0, -146100),
                               aa = c(0.999999570867702, 1.00000267419287),
                               ee = c(0.0167054504495442, 0.016854305837952),
                               inc = c(7.15495578299901, 7.14593294487681),
                               lph = c(27.323464691619, 26.1208300138006),
                               lan = c(179.999248773572, -179.58527791053),
                               arp = c(-152.675784081952, -154.293892075669),
                               mna = c(-2.45433424316846, 1.26740886389077)),
                          row.names = c(NA, -2L),
                          class = c("tbl_df", "tbl", "data.frame"))
  expect_snapshot(get_solution(astronomical_solution = ZB18a_head, quiet = TRUE))
})

test_that("get_solution() can load eccentricity solutions", {
  skip_on_cran() # because this needs to download the solutions

  pth <- withr::local_tempdir(pattern = "snvecR")
  withr::local_options(list(snvecR.cachedir = pth, width = 150))

  expect_snapshot(head(get_solution(astronomical_solution = "ZB20b",
                                    quiet = TRUE, force = TRUE)))

})

test_that("get_solution() can load full solutions", {
  skip_on_cran() # because this needs to download the solutions

  pth <- withr::local_tempdir(pattern = "snvecR")
  withr::local_options(list(snvecR.cachedir = pth, width = 150))

  expect_snapshot(head(get_solution(astronomical_solution = "full-ZB18a", quiet = FALSE, force = TRUE)),
                  # get rid of the cache directory printing in this snapshot because it differs between CIs and machines
                  transform = ~ gsub("^i The cache directory is '.*'.$",
                                     "i The cache directory is 'transformed-for-CI'.",
                                     .x))

  expect_snapshot(head(get_solution(astronomical_solution = "full-ZB18a", quiet = TRUE)))
})

test_that("get_solution() can load PT solutions", {
  skip_on_cran() # because this needs to download the solutions

  pth <- withr::local_tempdir(pattern = "snvecR")
  withr::local_options(list(snvecR.cachedir = pth, width = 150))

  expect_snapshot(head(get_solution(astronomical_solution = "PT-ZB18a(1,1)", quiet = FALSE, force = TRUE)),
                  # get rid of the cache directory printing in this snapshot because it differs between CIs and machines
                  transform = ~ gsub("^i The cache directory is '.*'.$",
                                     "i The cache directory is 'transformed-for-CI'.",
                                     .x))
  expect_snapshot(head(get_solution(astronomical_solution = "PT-ZB18a(1,1)", quiet = TRUE)))
})

test_that("get_solution() can load ZB23.Rxx solutions", {
  skip_on_cran() # because this needs to download the solutions

  pth <- withr::local_tempdir(pattern = "snvecR")
  withr::local_options(list(snvecR.cachedir = pth, width = 150))

  expect_snapshot(head(get_solution(astronomical_solution = "ZB23.R01", quiet = FALSE, force = TRUE)),
                  # get rid of the cache directory printing in this snapshot because it differs between CIs and machines
                  transform = ~ gsub("^i The cache directory is '.*'.$",
                                     "i The cache directory is 'transformed-for-CI'.",
                                     .x))
  expect_snapshot(head(get_solution(astronomical_solution = "ZB23.R01", quiet = TRUE)))
})
