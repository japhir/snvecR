test_that("get_solution() inputs are checked", {
  expect_error(get_solution(orbital_solution = "hoi"))
  expect_error(get_solution(orbital_solution = "La11"))
})

test_that("get_solution() works", {
  withr::local_options(width = 150)
  expect_snapshot(get_solution(orbital_solution = "ZB18a", quiet = TRUE) |> head())
  expect_snapshot(get_solution(orbital_solution = "ZB18a", quiet = FALSE, force = TRUE) |> head())
})
