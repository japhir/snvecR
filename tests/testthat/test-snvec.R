test_that("snvecR inputs are checked", {
  expect_error(snvec(orbital_solution = "hoi"))
  expect_error(snvec(orbital_solution = "La11"))
  expect_error(snvec(tend = -Inf))
  expect_error(snvec(tend = 5))
  expect_error(snvec(tres = -5))
  expect_error(snvec(tres = -Inf))
})

test_that("snvecR basic call works", {
  withr::local_options(width = 45)
  # I test a snapshot of the output
  # not the messages because they have timestamps
  expect_snapshot(snvec(quiet = TRUE) |>
                    dplyr::select(dplyr::all_of(c("age", "sx", "sy", "sx", "epl", "phi", "cp"))) |>
                    print(n = 100))
})
