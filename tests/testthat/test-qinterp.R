test_that("qinterp works", {
  expect_equal(
    qinterp(
      c(27.323464691619, 26.1208300138006, 24.6863325419433, 23.670354193149),
      # these are dat$lph[1:4]
      -146100, # that's dts
      -18262.5, # that's dx
      2 # that's m
    ),
    25.9415178298184
  )
})
