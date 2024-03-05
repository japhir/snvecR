test_that("get_solution() inputs are checked", {
  expect_error(get_solution(astronomical_solution = "hoi"))
  expect_error(get_solution(astronomical_solution = "full-La11"))
  expect_error(get_solution(astronomical_solution = data.frame(x = 1, y = 5))) # this runs via prepare_solution
})

test_that("get_solution() can load eccentricity solutions", {
  withr::local_options(width = 150)

  expect_snapshot(head(get_solution(astronomical_solution = "ZB20b",
                                    quiet = TRUE, force = TRUE)))
})

test_that("get_solution() can load full solutions", {
  withr::local_options(width = 150)

  expect_snapshot(head(get_solution(astronomical_solution = "full-ZB18a", quiet = FALSE, force = TRUE)),
                  # get rid of the cache directory printing in this snapshot because it differs between CIs and machines
                  transform = ~ gsub("^i The cache directory is '.*'.$",
                                     "i The cache directory is 'transformed-for-CI'.",
                                     .x))

  expect_snapshot(head(get_solution(astronomical_solution = "full-ZB18a", quiet = TRUE)))
})
