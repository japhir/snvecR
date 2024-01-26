test_that("unwrap works", {
  # the first 5 values of lan already have a flip
  # NOTE: I've also tried expect_identical to test machine precision, this is
  # NOT the same but only after 12 decimals
  expect_equal(unwrap(get_solution()$lan[1:6])[1:5],
               # output copied from C debugging section
               c(179.999248773572, 180.41472208947, 180.828937891654, 181.245975149062,
                 181.660845072088))
  # manually copied the first 5 values of orbitN-DE441 past run lph value
  lph <- c(-37.374247500332388, 105.957354230121595, 103.452633405997588,
           101.499240411608966, 99.018805879034474)
  expect_equal(unwrap(lph), lph)

  # this is lph in ZB18a at 472:475
  lph <- c(-179.038719967268, -179.7991103713, 179.458353781545, 178.80185739837)
  expect_equal(unwrap(lph),
               c(-179.038719967268, -179.7991103713, -180.541646218455, -181.19814260163))

})
