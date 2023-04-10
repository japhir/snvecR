test_that("snvecR inputs are checked", {
  expect_error(snvec(orbital_solution = "hoi"))
  expect_error(snvec(orbital_solution = "La11"))
  expect_error(snvec(tend = -Inf))
  expect_error(snvec(tend = 5))
  expect_error(snvec(tres = -5))
  expect_error(snvec(tres = -Inf))
})
