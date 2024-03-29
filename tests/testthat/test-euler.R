test_that("euler transform works", {
  expect_equal(
    as.vector(euler(c(1, 1, 1), 7.155, 75.5940, TRUE)),
    c(1.004640926031757, 0.074800407756640, 1.408936304004370)
  )
  expect_equal(
    as.vector(euler(c(-1, 0, 3), 7.155, 75.5940, TRUE)),
    c(-0.534164481136257, -2.447151567478811, 1.930315392074164)
  )
})
