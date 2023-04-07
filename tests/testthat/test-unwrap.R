test_that("unwrap works", {
  # the first 5 values of lan already have a flip
  expect_equal(unwrap(ZB18a$lan[1:6])[1:5],
               # output copied from C debugging section
               c(-179.585277910530010, -179.171062108345808,
                 -178.754024850938009, -178.339154927912205, -177.922959468805288))
})
