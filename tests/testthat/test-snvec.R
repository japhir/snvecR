test_that("snvecR inputs are checked", {
  expect_error(snvec(orbital_solution = "hoi"))
  expect_error(snvec(orbital_solution = "La11"))
  expect_error(snvec(tend = -Inf))
  expect_error(snvec(tend = 5))
  expect_error(snvec(tres = -5))
  expect_error(snvec(tres = -Inf))
  expect_error(snvec(output = "banaan"))
  expect_warning(snvec(ed = 1.5, tend = -1, tres = .5, quiet = TRUE))
  expect_warning(snvec(ed = 0.9, tend = -1, tres = .5, quiet = TRUE))
  expect_warning(snvec(td = 5, tend = -1, tres = .5, quiet = TRUE))
  expect_warning(snvec(td = -2, tend = -1, tres = .5, quiet = TRUE))
})

test_that("snvecR basic call works", {
  withr::local_options(width = 57)
  # I test a snapshot of the output
  expect_snapshot(snvec(tend = -49, # limit output to 50 kyr so it takes <5 s (CRAN check)
                        # specify default values explicitly in case they change in the future
                        ed = 1,
                        td = 0,
                        # return results at a very low resolution for speed
                        tres = 1,
                        # do not show info messages because some have timestamps
                        quiet = TRUE,
                        # provide output = all so that we can see if the vectors are doing well
                        output = "all") |>
                    # report only the columns of interest (the ones returned by the C-routine)
                    # this means that if I change my mind about which columns to report it doesn't matter,
                    # as long as the output of these columns remains the same.
                    dplyr::select(dplyr::all_of(c("age", "sx", "sy", "sz", "epl", "phi", "cp"))) |>
                    # print the full 100 rows to monitor changes
                    print(n = 50))
})
