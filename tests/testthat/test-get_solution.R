test_that("get_solution() inputs are checked", {
  expect_error(get_solution(orbital_solution = "hoi"))
  expect_error(get_solution(orbital_solution = "La11"))
  expect_error(get_solution(orbital_solution = data.frame(x = 1, y = 5))) # this runs via prepare_solution
})

test_that("get_solution() works", {
  withr::local_options(width = 150)
  expect_snapshot(get_solution(orbital_solution = "ZB18a", quiet = TRUE) |> head())
  expect_snapshot(get_solution(orbital_solution = "ZB18a", quiet = FALSE, force = TRUE) |> head(),
                  # get rid of the cache directory printing in this snapshot because it differs between CIs and machines
                  transform = ~ gsub("^i The cache directory is '.*'.$",
                                     "i The cache directory is 'transformed-for-CI'.",
                                     .x)
                  )
})
