#' Calculate Earth’s Obliquity and Precession in the Past
#'
#' `snvec()` computes climatic precession and obliquity (or tilt) from an
#' astronomical solution (AS) input and input values for dynamical ellipticity
#' (\eqn{E_{d}}{Ed}) and tidal dissipation (\eqn{T_{d}}{Td}). It solves a set
#' of ordinary differential equations.
#'
#' @param tend Final timestep in thousands of years (kyr).
#'   Defaults to `-1000` kyr.
#' @param ed Dynamical ellipticity \eqn{E_{d}}{Ed}, normalized to modern.
#'   Defaults to `1.0`.
#' @param td Tidal dissipation \eqn{T_{d}}{Td}, normalized to modern. Defaults
#'   to `0.0`.
# inherit astronomical_solution and quiet from get_solution
#' @inheritParams get_solution
#' @param os_ref_frame Character vector with the reference frame of the astronomical
#'   solution. Either `"HCI"` for heliocentric inertial reference frame or
#'   `"J2000"` for ecliptic J2000 reference frame. Defaults to `"HCI"` for
#'   `HNBody` output.
#' @param os_omt Longitude of ascending node of the solar equator relative to ECLIPJ2000.
#' @param os_inct Inclination of the solar equator relative to ECLIPJ2000.
# ODE solver parameters
#' @param tres Output timestep resolution in thousands of years (kyr). Defaults
#'   to `-0.4`. To determine the sign, think of from `0` to `tend` by timestep
#'   `tres`.
#' @param atol Numerical absolute tolerance passed to [deSolve::ode()]'s
#'   `atol`. Defaults to `1e-5`.
#' @param rtol Numerical relative tolerance passed to [deSolve::ode()]'s
#'   `rtol`. Defaults to `0`.
#' @param solver Character vector specifying the method passed to
#'   [deSolve::ode()]'s `method`. Defaults to `"vode"` for stiff problems with
#'   a variable timestep.
#' @param output Character vector with name of desired output. One of:
#'
#'   * `"nice"` (the default) A [tibble][tibble::tibble-package] with the
#'     columns `time`, `eei`, `epl`, `phi`, `cp`.
#'
#'   * `"full"` A [tibble][tibble::tibble-package] with all the computed and
#'     interpolated columns.
#'
#'   * `"ode"` A matrix with the output of the ODE solver.
#'
#' @details This is a re-implementation of the C-code in the supplementary
#'   information of Zeebe & Lourens (2022). The terms are explained in detail
#'   in Zeebe (2022).
#'
#' @section Reference Frames of Astronomical Solutions:
#'
#' NASA provides their asteroid and planet positions in the ecliptic J2000
#' reference frame, while long-term astronomical solution integrations are
#' often performed in the heliocentric inertial reference frame (HCI) or in the
#' inertial reference frame. This is to align the reference frame with the spin
#' vector of the Sun, making J2 corrections intuitive to implement.
#'
#' Obliquity is typically given in the ecliptic reference frame, so snvec
#' converts all outputs to J2000 if the `os_ref_frame` is equal to `"HCI"` and
#' does no transformations if it is already in `"J2000"`.
#'
#' For this, it uses \eqn{\Omega_{\odot} = 75.5940}{OMT = 75.5940} and
#' \eqn{i_{\odot} = 7.155}{INCT = 7.155} as in Zeebe (2017). You can overwrite
#' these defaults with `os_omt` and `os_inct` if desired.
#'
#' @section ODE Solver:
#'
#' Note that the different ODE solver algorithm we use (Soetaert et al.,
#' 2010) means that the R routine returns an evenly-spaced time grid, whereas
#' the C-routine has a variable time-step.
#' This means we need to explicitly set the stepsize `tres`.
#'
#' @returns `snvec()` returns different output depending on the `outputs` argument.
#'
#' If `output = "nice"` (the default), returns a
#' [tibble][tibble::tibble-package] with the following columns:
#'
# #'   * `t` Time \eqn{t} (days).
# #'
#'   * `time` Time in thousands of years (kyr).
#'
# I removed this earlier
# #'   * `eei` Astronomical solution's eccentricity \eqn{e}, interpolated to output
# #'   timescale (-).
# #'
#'   * `epl` Calculated Obliquity \eqn{\epsilon} (radians).
#'
#'   * `phi` Calculated Precession \eqn{\phi} (radians) from ECLIPJ2000.
#'
#'   * `cp` Calculated Climatic precession (-) as \eqn{e\sin\bar{\omega}}.
#'
#' where \eqn{\bar{\omega}} is the longitude of perihelion relative to the moving equinox.
#'
#' If `output = "all"` (for developers), additional columns are included,
#' typically interpolated to output timescale.
#'
#'   * `sx`, `sy`, `sz` The \eqn{x}{s}, \eqn{y}, and \eqn{z}-components of Earth's
#'   spin axis unit vector \eqn{\vec{s}}{s} in the heliocentric inertial
#'   reference frame.
#   this one is in HCI
#'
#'  See the source code for descriptions of all the intermediate computational
#'  steps.
#'
# #'
# #'   * `nnx`, `nny`, `nnz` The \eqn{x}, \eqn{y}, and \eqn{z}-components of the
# #'   unit normal vector \eqn{\vec{n}}{n}, normal to Earth's
# #'   instantaneous orbital plane.
# #   this one is in HCI
# #'
# #'   * `lphi` Unwrapped longitude of perihelion \eqn{\varpi} (radians).
# #'
# #'   * `lani` Unwrapped longitude of the ascending node \eqn{\Omega} (radians).
# #'
# #'   * `u` Spin axis unit vector \eqn{\vec{s}}{s} as a list-column.
# #'
# #'   * `nv` Unit normal vector to the orbital plane \eqn{\vec{n}}{n} as
# #'   a list-column.
# #'
# #'   * `up` Vector \eqn{\vec{u}'}{u'}, euler transform of
# #'   \eqn{\vec{u}}{u} to the instantaneous orbit plane (relative to
# #'   \eqn{\phi(t=0)=0} at J2000) as a list column.
# #   this one is in inertial frame ECLIPJ2000
#'
#' If `output = "ode"`, it will return the raw output of the ODE solver, which
#' is an object of class `deSolve` and `matrix`, with columns `time`, `sx`,
#' `sy`, and `sz`. This can be useful for i.e. [deSolve::diagnostics()].
#'
#' @seealso
#'
#' * [deSolve::ode()] from Soetaert et al., (2010) for the ODE solver that we
#'   use.
#'
#' * [get_solution()] Get astronomical solutions.
#'
#' @references
#'
#' Zeebe, R.E. (2017). Numerical Solutions for the Orbital Motion of the Solar
#'   System over the Past 100 Myr: Limits and New Results. _The Astronomical
#'   Journal_, 154(5), \doi{10.3847/1538-3881/aa8cce}.
#'
#' Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
#'  Paleocene–Eocene boundary age constrained by geology and astronomy.
#'  _Science_, 365(6456), 926–929. \doi{10.1126/science.aax0612}.
#'
#' Zeebe, R. E., & Lourens, L. J. (2022). A deep-time dating tool for
#'   paleo-applications utilizing obliquity and precession cycles: The role of
#'   dynamical ellipticity and tidal dissipation. _Paleoceanography and
#'   Paleoclimatology_, e2021PA004349. \doi{10.1029/2021PA004349}.
#'
#' Zeebe, R. E. (2022). Reduced Variations in Earth’s and Mars’ Orbital
#'   Inclination and Earth’s Obliquity from 58 to 48 Myr ago due to Solar
#'   System Chaos. _The Astronomical Journal_, 164(3),
#'   \doi{10.3847/1538-3881/ac80f8}.
#'
#' Wikipedia page on Orbital Elements:
#'   <https://en.wikipedia.org/wiki/Orbital_elements>
#'
#' Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer (2010). Solving
#'   Differential Equations in R: Package deSolve. Journal of Statistical
#'   Software, 33(9), 1–25. \doi{10.18637/jss.v033.i09}.
#'
#' @examples
#' \donttest{
#' \dontshow{
#' # set the cachedir to a temporary directory
#' pth <- withr::local_tempdir(pattern = "snvecR")
#' withr::local_options(snvecR.cachedir = pth)
#' }
#' # default call
#' snvec(tend = -1e3, ed = 1, td = 0, tres = -0.4)
#' }
#' @export
snvec <- function(tend = -1e3,
                  ed = 1,
                  td = 0,
                  astronomical_solution = "full-ZB18a",
                  os_ref_frame = "HCI",
                  os_omt = NULL, os_inct = NULL,
                  tres = -0.4,
                  atol = 1e-5,
                  rtol = 0,
                  solver = "vode",
                  quiet = FALSE,
                  output = "nice") {

  outputs <- c("nice", "all", "ode")
  if (!any(output == outputs)) {
    cli::cli_abort(c(
      "{.var output} must be one of {.or {.q {outputs}}}.",
      "x" = "You've supplied {.q {output}}."
    ))
  }

  ## tres
  if (sign(tres) != sign(tend)) {
    cli::cli_abort(c(
      "{.var tres} must be given in the same sign as {.var tend}.",
      "i" = "{.var tres} = {tres}",
      "i" = "{.var tend} = {tend}"
    ))
  }


  ## a quick dumb input test for now
  if (abs(tres) > abs(tend)) {
    cli::cli_abort(c(
      "abs({.var tres}) must be < abs({.var tend}).",
      "i" = "{.var tres} = {tres}",
      "i" = "{.var tend} = {tend}"
    ))
  }

  # this warning is too strict and kind of annoying
  ## if (ed < .998 | ed > 1.0005) {
  if (ed < .9 | ed > 1.1) {
    cli::cli_warn(c(
      "Dynamic ellipticity likely varied between 0.9980 and 1.0005 during the past 45 Ma!",
      "i" = "{.var ed} = {ed}",
      "*" = "See Zeebe & Lourens 2022 Pal&Pal <https://doi.org/10.1029/2021PA004349>."
    ))
  }

  if (td < 0 | td > 1.2) {
    cli::cli_warn(c(
      "Tidal dissipation likely varied between 0 and 1!",
      "i" = "{.var td} = {td}",
      "*" = "See Zeebe & Lourens 2022 Pal&Pal <https://doi.org/10.1029/2021PA004349>."
    ))
  }

  if (!any("data.frame" == class(astronomical_solution)) &&
        !grepl("^full-", astronomical_solution)) {
    cli::cli_abort(c(
      "Astronomical Solution must contain all orbital parameters",
      "i" = "Did you mean to specify {.q full-ZB18a}?"
    ))
  }

  hci_refs <- c("heliocentric intertial", "HCI")
  j2000_refs <- c("ecliptic", "ECLIPJ2000", "J2000")
  if (!any(os_ref_frame == c(hci_refs, j2000_refs))) {
    cli::cli_abort(c(
      "i" = "{.var os_ref_frame} must be one of {.or {.q {reference_frames}}}",
      "x" = "{.var os_ref_frame} = {.q {os_ref_frame}}"
    ))
  } else {
    # get rid of aliases
    if (any(os_ref_frame == hci_refs)) os_ref_frame <- "HCI"
    if (any(os_ref_frame == j2000_refs)) os_ref_frame <- "J2000"
  }

  # if the user specifies a value for os_omt and/or os_inct, they cannot
  # specify a non-default os_ref_frame
  if (os_ref_frame != "HCI") {
    if (!is.null(os_omt)) {
      cli::cli_abort(
        "Specified both {.var os_ref_frame} and {.var os_omt}."
      )
    }
    if (!is.null(os_inct)) {
      cli::cli_abort(
        "Specified both {.var os_ref_frame} and {.var os_inct}."
      )
    }
  }

  # if the reference frame is already J2000, just set the angles to 0
  # (this makes the codebase simpler, no need for ifelse statements throughout.)
  # but we will put in explicit ifelse's for optimisation
  if (os_ref_frame == "HCI" && is.null(os_omt) && is.null(os_inct)) {
    OMT <- 75.5940
    INCT <- 7.155
  } else if (os_ref_frame == "J2000") {
    if (is.null(os_omt)) OMT <- 0.
    if (is.null(os_inct)) INCT <- 0.
  }
  # or if the user has specified their own rotation, that's fine too!
  if (!is.null(os_omt)) OMT <- os_omt
  if (!is.null(os_inct)) INCT <- os_inct

  if (atol < 1e-12 | atol > 1e-3) {
    cli::cli_warn("Input absolute tolerance should be between 1e-12 and 1e-3.")
  }
  if (rtol > 1e-3) {
    cli::cli_warn("Input relative tolerance should be smaller than 1e-3.")
  }

  dat <- get_solution(astronomical_solution = astronomical_solution, quiet = quiet)

  if ((sign(tend) != sign(dat$time[2])) || (abs(tend) > max(abs(dat$time)))) {
    cli::cli_abort(c(
      "{.var tend} must fall within astronomical solution time.",
      "i" = "The astronomical solution {sign(dat$time[2])*max(abs(dat$time))}.",
      "x" = "{.var tend} = {tend}."
    ))
  }

  # message user about inputs
  if (!quiet) {
    startdate <- Sys.time()
    cli::cli_inform(c(
      "This is {VER}",
      "Richard E. Zeebe",
      "Ilja J. Kocken",
      "",
      "Integration parameters:",
      "*" = "{.var tend} = {.val {tend}} kyr",
      "*" = "{.var ed} = {.val {ed}}",
      "*" = "{.var td} = {.val {td}}",
      "*" = "{.var astronomical_solution} = {.val {if (any('data.frame' == class(astronomical_solution))) 'user provided' else astronomical_solution}}",
      "*" = "{.var os_ref_frame} = {.val {os_ref_frame}}",
      "*" = "{.var os_omt} = {if (is.null(os_omt)) 'defaulting to' else ''} {.val {OMT}}",
      "*" = "{.var os_inct} = {if (is.null(os_inct)) 'defaulting to' else ''} {.val {INCT}}",
      "*" = "{.var tres} = {.val {tres}} kyr",
      "*" = "{.var atol} = {.val {atol}}",
      "*" = "{.var rtol} = {.val {rtol}}",
      "*" = "{.var solver} = {.val {solver}}",
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
  dts <- dat$t[2] - dat$t[1] # difference in time in days

  ## initial values for the spin vector s
  ## [[file:snvec-3.7.5/snvec-3.7.5.c::=== finits()][finits()]]

  ## use finits to get initial conditions in transformed ECLIPJ2000
  omt <- OMT / R2D
  inct <- INCT / R2D
  ep0 <- EP0 / R2D
  cs <- cos(ep0)

  # first row of nn -> needs to be a vector
  # orbit normal at t=0
  ninit <- as.vector(as.matrix(dplyr::select(
    dplyr::filter(dat, .data$t == 0),
    tidyselect::all_of(c("nnx", "nny", "nnz")))))

  if (os_ref_frame != "J2000") {
    # transform n => n'
    np <- euler(ninit, inct, omt, TRUE)
  } else {
    np <- ninit
  }

  # solve quadratic equation for s0'y
  a <- np[2] * np[2] + np[3] * np[3]
  b <- -2 * cs * np[2]
  c <- cs * cs - np[3] * np[3]

  # initial position of spin vector in Earth's mean equator frame at t=0 (J2000) = [0, 0, 1]
  s0p <- c(NA, NA, NA)
  s0p[2] <- (-b + sqrt(b * b - 4 * a * c)) / (2 * a)
  s0p[3] <- sqrt(1 - s0p[2] * s0p[2])
  s0p[1] <- 0
  as.matrix(s0p)

  if (os_ref_frame != "J2000") {
    # transform s0' to s0
    s0 <- euler(s0p, inct, omt, FALSE)
  } else {
    s0 <- s0p
  }

  ## set the deSolve state
  state <- c(
    sx = s0[1],
    sy = s0[2],
    sz = s0[3]
  )

  ## define deSolve parameters
  parameters <- c(
    k0d = k0d,
    wdw = wdw,
    ndn = ndn
  )

  ## the differential equations
  ## [[file:snvec-3.7.5/snvec-3.7.5.c::=== derivs() ][derivs()]]
  # derivatives. RHS of DEQs for spin vector s = y
  eqns <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      # the the index of the current timestep in the astronomical solution
      m <- min(round(abs(t / dts) + 1), nrow(dat))
      # quickly interpolate astronomical solution input value at the current
      # timestep
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
      # they're also equally fast!
      # we keep it the same as the C-routine for now.

      ff <- 1 / sqrt(ff * ff * ff)
      kb <- k0d * (1 + 1 * wdw * t) * (ff + BET0 * (1 + 2 * ndn * t))

      fac <- FGP * kb * (ddi * (ppi * sx - qqi * sy) + cci * sz)

      dX <- fac * (cci * sy + ddi * qqi * sz)
      dY <- fac * (-cci * sx + ddi * ppi * sz)
      dZ <- -fac * (qqi * sx + ppi * sy) * ddi

      # EPSDOT: tidal effect on obliquity
      # if (epsdot) {
      # orbit normal at ti = m is nn[j][m]
      # dot product and dydt:
      # dotab = cos(epl), sqrt = sin(epl=acos(dotab))

      ## dotab = s[1]*nn[1][m]+s[2]*nn[2][m]+s[3]*nn[3][m];
      ## tmp = tdg*EPSDOT*D2S/sqrt(1.-dotab*dotab);
      ## yp[1] += tmp*(nn[1][m] - dotab*s[1]);
      ## yp[2] += tmp*(nn[2][m] - dotab*s[2]);
      ## yp[3] += tmp*(nn[3][m] - dotab*s[3]);
      # }
      list(c(dX, dY, dZ))
    }) # end 'with(as.list( ...
  }

  ## a linear sequence of steps
  times <- seq(0, tend * KY2D,
               by = tres * KY2D
               # ~ tres = 0.4 is the average diff in the C-output
               # snv_sout$time |> diff() |> median() = 0.396013
               )

  ## solve it
  ## [[file:snvec-3.7.5/snvec-3.7.5.c::%%% solver][odeint()]]
  out <- deSolve::ode(
    y = state,
    times = times,
    func = eqns,
    parms = parameters,
    method = solver,
    atol = atol,
    rtol = rtol
    # perhaps I should pass ... here?
  )

  ## print the final values for s
  ## This is at t = tend, it's going back from 0 to -time by default.
  fin <- out[nrow(out), ]
  u <- as.vector(c(fin[2], fin[3], fin[4]))

  if (anyNA(fin)) {
    cli::cli_abort(c("The procedure returned NA for spin vector s!",
                     "i" = "{u}"))
  }

  if (!quiet) {
    cli::cli_inform(c(
      "Final values:",
      "*" = "s[1][2][3]: {.val {c(fin[2], fin[3], fin[4])}}",
      ## "*" = "s-error = |s|-1: {sqrt(pracma::dot(u, u))-1}"
      "*" = "s-error = |s|-1: {.val {sqrt(u[1]^2+u[2]^2+u[3]^2)-1}}"
    ))
  }

  # return ODE output if desired
  if (output == "ode") {
    return(out)
  }

  ## interpolate the full astronomical solution onto output timescale
  fin <- dplyr::mutate(
    dplyr::rename(t = "time",
                  tibble::as_tibble(out)),
    time = .data$t / KY2D, .after = "t")
  fin <- dplyr::mutate(fin,
                       nnx = approxdat(dat, "nnx")(.data$t),
                       nny = approxdat(dat, "nny")(.data$t),
                       nnz = approxdat(dat, "nnz")(.data$t),
                       eei = approxdat(dat, "ee")(.data$t),
                       inci = approxdat(dat, "inc")(.data$t),
                       lphi = approxdat(dat, "lphu")(.data$t),
                       lani = approxdat(dat, "lanu")(.data$t)
                       )

  fin <- dplyr::mutate(fin,
      # calculate obliquity
      ## tmp = pracma::dot(.data$u, .data$nv), # get rid of dependency
      tmp = .data$sx * .data$nnx + .data$sy * .data$nny + .data$sz * .data$nnz,
      # calculate obliquity
      epl = acos(.data$tmp))

  ## fin <- dplyr::mutate(dplyr::rowwise(fin),
    # for each row, NOTE this makes it very slow!!
  ## fin <- dplyr::mutate(fin,
    # extract sx, sy, sz, and nnx, nny, nnz
      # create list columns of vectors
      ## u = list(matrix(c(.data$sx, .data$sy, .data$sz), ncol = 1, nrow = 3)),
      ## nv = list(matrix(c(.data$nnx, .data$nny, .data$nnz), ncol = 1, nrow = 3)),
                       ## u = purrr::pmap(list(.data$sx, .data$sy, .data$sz), c))
                       ## u = list(c(.data$sx, .data$sy, .data$sz)))
  ## fin <- dplyr::mutate(dplyr::ungroup(fin), # end rowwise

  fin$u <- Map(c, fin$sx, fin$sy, fin$sz)
  # Map is slightly (but not a lot) faster than purrr::map
  fin$up <- Map(euler, s = fin$u, inc = fin$inci / R2D, lan = fin$lani / R2D)
  # coords: relative to phi(t=0)=0 at J2000
  fin$up <- Map(euler, s = fin$up, inc = 0, lan = -(fin$lani + OMT) / R2D + pi / 2)

  fin <- dplyr::mutate(fin,
      ## calculate axial precession
      # get 2nd and 1st value of up
      phi = purrr::map_dbl(.data$up, .f = \(x) atan2(x[2], x[1])),
      # normalize to first value of phi
      phi = .data$phi - dplyr::first(.data$phi),
      # calculate climatic precession
      cp = .data$eei * sin((.data$lphi + OMT) / R2D - .data$phi)
    )

  ## message user about final values
  if (!quiet) {
    stopdate <- Sys.time()
    cli::cli_inform(
      c("Final values:",
        "*" = "obliquity: {.val {fin[nrow(fin), 'epl']}} rad",
        "*" = "precession: {.val {fin[nrow(fin), 'phi']}} rad",
        "i" = "stopped at {.q {stopdate}}",
        "i" = "total duration: {.val {round(stopdate - startdate, 2)}}"
        )
    )
  }

  # final cleanup
  # we transform the deSolve parameters into simple numeric columns
  # this is so they work better with things like bind_rows etc. via vctrs
  fin <- dplyr::mutate(fin,
                       dplyr::across(
                         tidyselect::all_of(c("t", "time", "sx", "sy", "sz", "epl")),
                         as.numeric))

  if (output == "all") {
    return(fin)
  }

  if (output == "nice") {
    return(dplyr::select(fin,
                         tidyselect::all_of(c("time", "epl", "phi", "cp"))))
  }
}
