test_that("get_solution() inputs are checked", {
  expect_error(get_solution(orbital_solution = "hoi"))
  expect_error(get_solution(orbital_solution = "La11"))
})

test_that("get_solution() works", {
  expect_snapshot(get_solution(orbital_solution = "ZB18a") |> head())
})
