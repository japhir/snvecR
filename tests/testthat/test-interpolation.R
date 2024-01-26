test_that("approxdat() works", {
  # create very simple test dataframe
  dat <- tibble(t = seq(0, -10, -.2),
                y = sin(seq(-2 * pi, 2 * pi, length.out = length(t))))
  # test the vectorized method
  # that does interpolation for the whole dataset at once
  # with the target time exactly half-way between rows 5 and 6
  expect_equal(approxdat(dat, "y")((dat$t[6] + dat$t[5]) / 2),
               (dat$y[6] + dat$y[5]) / 2)
  ## dat |>
  ##   ggplot(aes(x = t, y = y)) +
  ##   geom_line() + geom_point() +
  ##   annotate("point",
  ##            x = dat$t + .1,
  ##            y = approxdat(dat, "y")(dat$t + .1),
  ##            size = 5, colour = "red")
  fut <- tibble(t = seq(0, 10, .2),
                y = seq(0, 1, length.out = length(t)))
  expect_equal(approxdat(fut, "y")((fut$t[6] + fut$t[5]) / 2),
               (fut$y[6] + fut$y[5]) / 2)
})

test_that("qinterp() works", {

  # "normal" example where we go from 0 to -t
  dat <- tibble(t = seq(0, -10, -2),
                y = sin(seq(.3 * pi, pi, length.out = length(t))))
  # do the same interpolation but using qinterp
  dts <- dat$t[2] - dat$t[1] # diff between timesteps, -0.2
  # target timestep half-way between rows 5 and 6
  tmv <- (dat$t[6] + dat$t[5]) / 2 # -0.9
  # index integer m, this is 5
  m <- min(round(abs(tmv / dts)) + 1, nrow(dat))
  ## dx <- dat$t[m] - t # i've experimented with this as well
  dx <- tmv - dat$t[m]
  ## # visualize the fake data
  ## dat |>
  ##   ggplot(aes(x = t, y = y)) +
  ##   geom_line() + geom_point() +
  ##   annotate("point", x = c(tmv, tmv),
  ##            y = c((dat$y[6] + dat$y[5]) / 2,
  ##                  qinterp(dat$y, dts, dx, m)),
  ##            size = c(5, 3), colour = c("darkgreen", "orange")) +
  ##   geom_rug(aes(x = t, y = NULL)) +
  ##   geom_rug(data = tibble(t = tmv, y = NA_real_), colour = "red", linewidth = 5)
  expect_equal(qinterp(dat$y, dts, dx, m),
               (dat$y[6] + dat$y[5]) / 2)

  # "future" example where we go from 0 to +t
  fut <- tibble(t = seq(0, 10, 2),
                y = sin(seq(.3 * pi, pi, length.out = length(t))))
  dts <- fut$t[2] - fut$t[1] # diff between timesteps
  tmv <- (fut$t[6] + fut$t[5]) / 2
  m <- min(round(abs(tmv / dts)) + 1, nrow(fut))
  dx <- tmv - fut$t[m]
  ## # visualize the fake data
  ## fut |>
  ##   ggplot(aes(x = t, y = y)) +
  ##   geom_line() + geom_point() +
  ##   annotate("point", x = c(tmv, tmv),
  ##            y = c((fut$y[6] + fut$y[5]) / 2, qinterp(fut$y, dts, dx, m)),
  ##            size = c(5, 3), colour = c("darkgreen", "orange")) +
  ##   geom_rug(aes(x = t, y = NULL)) +
  ##   geom_rug(data = tibble(t = tmv, y = NA_real_), colour = "red", linewidth = 5)

  expect_equal(qinterp(fut$y, dts, dx, m),
               (fut$y[6] + fut$y[5]) / 2)

  # more specific tests to see if it does the same as the C code
  # the first 4 values of ee
  ts <- c(0.000000000000000, -146100.000000000000000, -292200.000000000000000, -438300.000000000000000)
  ee <- c(0.016705450449544, 0.016854305837952, 0.017066750754228, 0.017191712714982)
  # desired timesteps
  tmv <- c(0.000000000000000, -137939.394901296764147, -276729.110959491459653, -414382.313519898743834)
  # change in timestep in input astronomical solution
  dts <- ts[2] - ts[1]
  # diff between desired time and os time
  ## dx <- tmv[2] - ts[2]
  dx <- 8160.605098703235853
  m <- 2
  mm <- 1
  ## dy <- ee[mm] - ee[m]
  dy <- -0.000148855388408
  ## yi <- ee[m] + dy * abs(dx) / abs(ds)
  # target
  yi <- 0.016845991327058
  ## tibble(t = ts, y = ee) |>
  ##   ggplot(aes(x = t, y = y)) +
  ##   geom_point() +
  ##   geom_line() +
  ##   annotate("point", x = c(tmv[2], tmv[2]),
  ##            y = c(yi, qinterp(ee, dts, dx, m)),
  ##            size = c(5, 3), colour = c("darkgreen", "orange")) +
  ##   geom_rug(data = tibble(t = ts, y = NA_real_)) +
  ##   geom_rug(data = tibble(t = tmv, y = NA_real_), colour = "red", linewidth = 5)

  expect_equal(qinterp(y = ee, ds = dts, dx = dx, m = m), yi)
  # does it still work if we flip time?
  expect_equal(qinterp(y = ee, ds = -dts, dx = -dx, m = m), yi)

})

test_that("qinterp and approxdat give the same results", {
  dat <- tibble(t = seq(0, -10, -.2),
                y = sin(seq(-2 * pi, 2 * pi, length.out = length(t))))
  expect_equal(approxdat(dat, "y")((dat$t[6] + dat$t[5]) / 2),
               qinterp(y = dat$y, ds = dat$t[2] - dat$t[1],
                       dx = ((dat$t[6] + dat$t[5]) / 2) - dat$t[5],
                       m = 5))
})
