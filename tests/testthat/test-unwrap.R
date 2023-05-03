test_that("unwrap works", {
  # the first 5 values of lan already have a flip
  # NOTE: I've also tried expect_identical to test machine precision, this is
  # NOT the same but only after 12 decimals
  expect_equal(unwrap(get_ZB18a()$lan[1:6])[1:5],
               # output copied from C debugging section
               c(179.999248773572, 180.41472208947, 180.828937891654, 181.245975149062,
                 181.660845072088))
})
