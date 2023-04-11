#' Calculate precession and obliquity from an OS, Td, and Ed.
#'
#' Computes climatic precession and obliquity (or tilt) from an orbital
#' solution (OS) input and input values for dynamical ellipticity (Ed) and
#' tidal dissipation (Td). It fits a set of ordinary differential equations.
#'
#' @param tend The final timestep in -kyr. Defaults to `-1000` kyr.
#' @param ed Dynamical ellipticity, normalized to modern. Defaults to `1.0`.
#' @param td Tidal dissipation, normalized to modern. Defaults to `0.0`.
#' @param orbital_solution Character vector with the name of the orbital
#'   solution to use. One of `"ZB18a"` (default) or `"La11"` (not yet
#'   implemented!).
#' @param tres The output timestep resolution in kyr. Defaults to `0.4`.
#' @param tolerance The numerical tolerance passed to [deSolve::ode()]'s `rtol`
#'   and `atol` arguments. Defaults to 1e-7.
#' @param quiet Be quiet?
#'   If `TRUE`, hide info messages.
#'   If `FALSE` (the default) print info messages and timing.
#' @returns A [tibble][tibble::tibble-package] with all the computed results (see Details).
#'
#' @author Ilja J. Kocken and Richard E. Zeebe
#'
#' @details
#' The output is a [tibble][tibble::tibble-package] with the following columns:
#'   * `time` Time in years.
#'   * `sx`, `sy`, `sz` Input vector s.
#'   * `age` Age in thousands of years ago (ka).
#'   * `nnx`, `nny`, `nnz` The euler transform of the input vector s.
#'   * `eei`, `lphi`, `lani` The orbital solution's eccentricity, unwrapped
#'      long periapse, and unwrapped long ascending node, interpolated to
#'      resulting timesteps.
#'   * `tmp` A temporary value that we used for the calculation.
#'   * `epl` The `acos(tmp)`.
#'   * `u` The input vector s as a list-column.
#'   * `nv` The vector n as a list-column.
#'   * `up` Vector u', with coordinates relative to phi(t=0) at J2000
#'   * `phi` The true anomaly.
#'   * `cp` Climatic precession.
#'
#' @seealso
#'   [deSolve::ode()] for the ODE solver that we use.
#'
#' @references
#'
#' Zeebe, R. E. (2022). Reduced Variations in Earth’s and Mars’
#'   Orbital Inclination and Earth’s Obliquity from 58 to 48 Myr ago due to
#'   Solar System Chaos. _The Astronomical Journal_, 164(3), 107.
#'   <https://doi.org/10.3847/1538-3881/ac80f8>
#'
#' Wikipedia page on Orbital Elements: <https://en.wikipedia.org/wiki/Orbital_elements>
#'
#' @examples
#' # default call
#' snvec()
#'
#' # a quick one with few timesteps, low resolution, high tolerance
#' snvec(-1e2, 1, 0, orbital_solution = "ZB18a", tres = 1, tolerance = 1e-4)
#' @export
snvec <- function(tend = -1e3,
                  ed = 1,
                  td = 0,
                  orbital_solution = "ZB18a",
                  tres = 0.4,
                  tolerance = 1e-7,
                  quiet = FALSE) {
  ## select the desired orbital solution
  solutions <- c("ZB18a", "La11")
  if (!orbital_solution %in% solutions) {
    cli::cli_abort(c("{.var orbital_solution} must be one of: {.or {.q {solutions}}}",
      "x" = "You've supplied {.q {orbital_solution}}"
    ))
  }
  if (orbital_solution == "ZB18a") {
    dat <- snvecR::ZB18a
  }
  if (orbital_solution == "La11") {
    ## dat <- snvecR::La11
    cli::cli_abort(c("Orbital solution: La11 currently not supported.",
      "i" = "Pull requests welcome."
    ))
  }

  ## tend
  if (tend >= 0) {
    cli::cli_abort(c("{.var tend} must be < 0",
      "x" = "You've supplied {tend}"
    ))
  }
  if (tend < min(dat$t / KY2D)) {
    cli::cli_abort(c("{.var tend} must be > the orbital solution {min(dat$t)/KY2D}",
      "x" = "You've supplied {tend}"
    ))
  }

  ## tres
  ## a quick dumb input test for now
  if (tres < tend) {
    cli::cli_abort(c("{.var tres} must be < {.var tend}.",
      "i" = "{.var tres} = {tres}",
      "i" = "{.var tend} = {tend}"
    ))
  }

  if (tres <= 0) {
    cli::cli_abort(c("{.var tres} must be > 0.",
                     "i" = "{.var tres} = {tres}"))
  }

  if (ed < .998 | ed > 1.0005) {
    cli::cli_warn(c("!" = "Dynamic ellipticity likely varied between 0.9980 and 1.0005 during the past 45 Ma!",
                    "i" = "{.var ed} = {ed}",
                    "*" = "See Zeebe & Lourens 2022 Pal&Pal <https://doi.org/10.1029/2021PA004349>"))
  }

  if (td < 0 | td > 1) {
    cli::cli_warn(c("!" = "Tidal dissipation likely varied between 0 and 1!",
                    "i" = "{.var td} = {td}",
                    "*" = "See Zeebe & Lourens 2022 Pal&Pal <https://doi.org/10.1029/2021PA004349>"))
  }

  if (tolerance > 1e-3 | tolerance < 1e-12) {
    cli::cli_warn("Input tolerance should be between 1e-3 and 1e-12.")
  }

  # message user about inputs
  if (!quiet) {
    startdate <- lubridate::now()
    cli::cli_inform(c(
      ## "This is {VER}",
      ## "Richard E. Zeebe",
      ## "Ilja J. Kocken",
      "Integration parameters:",
      "*" = "{.var tend} = {tend} kyr",
      "*" = "{.var ed} = {ed}",
      "*" = "{.var td} = {td}",
      "*" = "{.var orbital_solution} = {.q {orbital_solution}}",
      "*" = "{.var tres} = {tres} kyr",
      "*" = "{.var tolerance} = {tolerance}",
      "i" = "started at {.q {startdate}}"
    ))
  }

  ## calculate global vars ndn, wdw, k0d from Td and Ed
  ## [[file:snvec-3.7.5/snvec-3.7.5.c::=== fedtd() ][fedtd()]]

  # as a function of ed, td
  k0d <- ((3 / 2) * GM * ED0 * ed / (OM * AU3)) * D2S # 1/s => 1/d
  k0b0 <- k0d * (1 + BET0)
  ndn <- -4.6e-18 * D2S * td # 1/s => 1/d
  wdw <- 51 * ndn * NW0 # Lambeck80, see PTman
  tdg <- td # global Td
  dts <- dat$t[2] - dat$t[1] # difference in time

  ## initial values for the spin vector s
  ## [[file:snvec-3.7.5/snvec-3.7.5.c::=== finits()][finits()]]

  ## use finits to get initial conditions in transformed ECLIPJ2000
  omt <- OMT / R2D
  inct <- INCT / R2D
  ep0 <- EP0 / R2D
  cs <- cos(ep0)

  # first row of nn -> needs to be a vector
  # orbit normal at t=0
  ninit <- dat |>
    dplyr::filter(.data$t == 0) |>
    dplyr::select(tidyselect::all_of(c("nnx", "nny", "nnz"))) |>
    as.matrix() |>
    as.vector()

  # transform n => n'
  np <- euler(ninit, inct, omt, TRUE)

  # solve quadratic equation for s0'y
  a <- np[2] * np[2] + np[3] * np[3]
  b <- -2 * cs * np[2]
  c <- cs * cs - np[3] * np[3]

  s0p <- c(NA, NA, NA)
  s0p[2] <- (-b + sqrt(b * b - 4 * a * c)) / (2 * a)
  s0p[3] <- sqrt(1 - s0p[2] * s0p[2])
  s0p[1] <- 0
  as.matrix(s0p)

  # transform s0' to s0
  s0 <- euler(s0p, inct, omt, 0)

  ## set the deSolve state
  state <- c(
    sx = s0[1],
    sy = s0[2],
    sz = s0[3]
  )

  ## define deSolve parameters
  parameters <- c(
    ed = ed,
    td = td,
    k0d = k0d,
    wdw = wdw,
    ndn = ndn
  )

  ## the differential equations
  ## [[file:snvec-3.7.5/snvec-3.7.5.c::=== derivs() ][derivs()]]
  # derivatives. RHS of DEQs for spin vector s = y
  eqns <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      m <- min(round(abs(t / dts) + 1), nrow(dat))

      dx <- t - dat$t[m]
      qqi <- qinterp(dat$qq, dts, dx, m)
      ppi <- qinterp(dat$pp, dts, dx, m)
      cci <- qinterp(dat$cc, dts, dx, m)
      ddi <- qinterp(dat$dd, dts, dx, m)

      # 1/(1-e^2)^3/2 term
      # add interpolation
      ## hhi <- qinterp(dat$hh,dts,dx,m)
      ## kki <- qinterp(dat$kk,dts,dx,m)
      ## ff <- (1 - hhi * hhi - kki * kki)

      # shouldn't I also interpolate hh and kk? -> see above
      ff <- (1 - dat$hh[m] * dat$hh[m] - dat$kk[m] * dat$kk[m])
      # i've tried both, gives identical results if I use the prescribed timesteps.
      # they're also equally fast! so let's go with my own which I think is better.
      # it might be the cause of numerical diffs between C and R? try without again

      ff <- 1 / sqrt(ff * ff * ff)
      kb <- k0d * (1 + 1 * wdw * t) * (ff + BET0 * (1 + 2 * ndn * t))

      fac <- FGP * kb * (ddi * (ppi * sx - qqi * sy) + cci * sz)

      dX <- fac * (cci * sy + ddi * qqi * sz)
      dY <- fac * (-cci * sx + ddi * ppi * sz)
      dZ <- -fac * (qqi * sx + ppi * sy) * ddi

      # EPSDOT
      ## dotab = s[1]*nn[1][m]+s[2]*nn[2][m]+s[3]*nn[3][m];
      ## tmp = tdg*EPSDOT*D2S/sqrt(1.-dotab*dotab);
      ## yp[1] += tmp*(nn[1][m] - dotab*s[1]);
      ## yp[2] += tmp*(nn[2][m] - dotab*s[2]);
      ## yp[3] += tmp*(nn[3][m] - dotab*s[3]);

      list(c(dX, dY, dZ))
    }) # end 'with(as.list( ...
  }

  ## a linear sequence of steps

  times <- seq(0, tend * KY2D,
               ## length.out = 2523L # the length of the C-output
               by = -tres * KY2D
               # ~ tres = 0.4 is the average diff in the C-output
               # snv_sout$time |> diff() |> median() = 0.396013
  )

  ## solve it
  ## [[file:snvec-3.7.5/snvec-3.7.5.c::%%% solver][odeint()]]
  ## print(system.time(
    ## microbenchmark::microbenchmark(
    out <- deSolve::ode(
      y = state,
      times = times,
      func = eqns,
      parms = parameters,
      method =
      ## "lsoda"# = default, chooses stiff/nonstiff automatically starting non-stiff
      # "ode23" # = non-stiff, variable time-step
      ## "ode45" # = stiff, variable time-step
      # radau #= stiff/non-stiff
        "bdf", # = stiff
      ## "daspk", # = very stiff
      # play around with machine precision: default is 1e-6
      ## rtol = 1e-5, atol = 1e-5 # rougher = faster?
      rtol = tolerance, atol = tolerance # based on EPSLVR = 1e-7
      ## rtol = 1e-12, atol = 1e-12
    )
  ## ))
  ## )

  ## print the final values for s
  ## This is at t = tend, it's going back from 0 to -time
  fin <- out[nrow(out), ]
  u <- as.vector(c(fin[2], fin[3], fin[4]))

  if (!quiet) {
    cli::cli_inform(c(
      "Final values:",
      "*" = "s[1][2][3]: {paste(fin[2], fin[3], fin[4])}",
      "*" = "s-error = |s|-1: {sqrt(abs(pracma::dot(u, u)))-1}"
    ))
  }

  ## interpolate the orbital solution
  ## back onto output timescale
  fin <- out |>
    tibble::as_tibble() |>
    dplyr::mutate(
      age = -.data$time / KY2D,
      nnx = approxdat(dat, "nnx")(.data$time),
      nny = approxdat(dat, "nny")(.data$time),
      nnz = approxdat(dat, "nnz")(.data$time),
      eei = approxdat(dat, "ee")(.data$time),
      inci = approxdat(dat, "inc")(.data$time),
      lphi = approxdat(dat, "lphu")(.data$time),
      lani = approxdat(dat, "lanu")(.data$time)
    )

  ## calculate obliquity
  fin <- fin |>
    # calculate the dotproduct, Richard's vvdot
    dplyr::mutate(
      tmp = .data$sx * .data$nnx + .data$sy * .data$nny + .data$sz * .data$nnz,
      epl = acos(.data$tmp)
    )

  ## calculate precession and climatic precession
  fin <- fin |>
    # for each row, NOTE this makes it very slow!!
    dplyr::rowwise() |>
    # create list columns of matrices
    # extract sx, sy, sz, and nnx, nny, nnz
    dplyr::mutate(
      u = list(matrix(c(.data$sx, .data$sy, .data$sz), ncol = 1, nrow = 3)),
      nv = list(matrix(c(.data$nnx, .data$nny, .data$nnz), ncol = 1, nrow = 3)),
      # coords: fixed HCI => moving orbit plane
      up = list(euler(.data$u, .data$inci / R2D, .data$lani / R2D, 0)),
      # coords: relative to phi(t=0)=0 at J2000
      up = list(euler(.data$up, 0, -(.data$lani + OMT) / R2D - pi / 2, 0)),
      # get 2nd and 1st column of up
      phi = map2_dbl(.data$up[2, ], .data$up[1, ], atan2)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      # normalize to first value of phi
      tmp = first(.data$phi),
      phi = .data$phi - .data$tmp,
      cp = .data$eei * sin((.data$lphi + OMT) / R2D - .data$phi),
      # phi is now mapped between 0 and 2 pi, whereas RZ's output is wrapped
      # between -pi and pi
      phi = ifelse(.data$phi > pi, .data$phi - 2 * pi, .data$phi)
    )

  ## message user about final values
  if (!quiet) {
    cli::cli_inform(
      c("Final values:",
        "*" = "obliquity: {fin[nrow(fin), 'epl']} rad",
        "*" = "precession: {fin[nrow(fin), 'phi']} rad",
        "i" = "stopped at {.q {lubridate::now()}}",
        "i" = "total duration:  {lubridate::as.duration(round(lubridate::now() - startdate, 2))}"
      )
    )
  }

  # return fin
  fin |>
    # we transform the deSolve parameters into simple numeric columns
    # this is so they work better with things like bind_rows etc. via vctrs
    dplyr::mutate(dplyr::across(
      all_of(c("time", "sx", "sy", "sz", "age", "epl")),
      as.numeric))
}
