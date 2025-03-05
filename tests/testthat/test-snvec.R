test_that("snvec() invalid inputs throw errors", {
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

  # we test the astronomical_solution inputs for the get_solution helper in stead
  expect_error(snvec(tend = -Inf, quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(tres = -Inf, quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(tend = -1000, tres = 0.4, quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(tend = 1000, tres = -0.4, quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(tend = 1000, tres = 0.4, quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(tend = -1000, tres = -0.4, quiet = TRUE,
                     astronomical_solution = dplyr::mutate(ZB18a_head, t = -t)))

  expect_error(snvec(output = "banaan", quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(solver = "banaan", quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  # same as get_solution
  expect_error(snvec(astronomical_solution = "hoi", quiet = TRUE))
  expect_error(snvec(os_ref_frame = "hoi", quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(os_ref_frame = "J2000", os_omt = 5, quiet = TRUE,
                     astronomical_solution = ZB18a_head))
  expect_error(snvec(os_ref_frame = "J2000", os_inct = 7, quiet = TRUE,
                     astronomical_solution = ZB18a_head))

})

test_that("snvec() inputs outside of likely range throw warnings", {
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

  # these are a bit annoying
  expect_warning(snvec(ed = 1.1001,
                       tend = -.4, tres = -.4, quiet = TRUE,
                       astronomical_solution = ZB18a_head))
  expect_warning(snvec(ed = 0.8999,
                       tend = -.4, tres = -.4, quiet = TRUE,
                       astronomical_solution = ZB18a_head))
  expect_warning(snvec(td = 1.2001,
                       tend = -.4, tres = -.4, quiet = TRUE,
                       astronomical_solution = ZB18a_head))
  expect_warning(snvec(td = -0.0001,
                       tend = -.4, tres = -.4, quiet = TRUE,
                       astronomical_solution = ZB18a_head))
  expect_warning(snvec(atol = 0.9e-12,
                       tend = -.4, tres = -.4, quiet = TRUE,
                       astronomical_solution = ZB18a_head))
  expect_warning(snvec(atol = 1.1e-3,
                       tend = -.4, tres = -.4, quiet = TRUE,
                       astronomical_solution = ZB18a_head))
  expect_warning(snvec(rtol = 1.1e-3,
                       tend = -.4, tres = -.4, quiet = TRUE,
                       astronomical_solution = ZB18a_head))
})

test_that("snvec() output columns are correct", {
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

  # we have the desired columns
  expect_equal(colnames(
    snvec(tend = -.4, tres = -.4, quiet = TRUE, output = "nice",
          astronomical_solution = ZB18a_head)),
    c("time", "epl", "phi", "lpx", "cp"))

  # we have a deSolve matrix
  expect_equal(
    class(snvec(tend = -.4, tres = -.4, quiet = TRUE, output = "ode",
                astronomical_solution = ZB18a_head)),
    c("deSolve", "matrix"))
})

test_that("snvec() works", {
  skip_on_cran()

  pth <- withr::local_tempdir(pattern = "snvecR")
  withr::local_options(list(snvecR.cachedir = pth, width = 100))

  # I test a snapshot of the output
  expect_snapshot(
    # print the full 100 rows to monitor changes
    print(
      dplyr::select(
        snvec(tend = -49, # limit output to 50 ka so it takes <5 s (CRAN check)
              # specify default values explicitly in case they change in the future
              ed = 1,
              td = 0,
              # return results at a very low resolution for speed
              tres = -1,
              os_ref_frame = "HCI",
              # do not show info messages because some have timestamps
              quiet = TRUE,
              # provide output = all so that we can see if the vectors are doing well
              output = "all"),
        # report only the columns of interest (the ones returned by the C-routine)
        # this means that if I change my mind about which columns to report it doesn't matter,
        # as long as the output of these columns remains the same.
        tidyselect::all_of(c("time", "sx", "sy", "sz", "epl", "phi", "lpx", "cp"))),
      n = 50))

    expect_snapshot(
    # print the full 100 rows to monitor changes
    print(
      dplyr::select(
        snvec(tend = -49, # limit output to 50 ka so it takes <5 s (CRAN check)
              # specify default values explicitly in case they change in the future
              ed = 1,
              td = 0,
              # return results at a very low resolution for speed
              tres = -1,
              # do not show info messages because some have timestamps
              quiet = TRUE,
              # provide output = all so that we can see if the vectors are doing well
              os_ref_frame = "J2000",
              # set os_ref_frame to J2000
              output = "all"),
        # report only the columns of interest (the ones returned by the C-routine)
        # this means that if I change my mind about which columns to report it doesn't matter,
        # as long as the output of these columns remains the same.
        tidyselect::all_of(c("time", "sx", "sy", "sz", "epl", "phi", "lpx", "cp"))),
      n = 50))

})
