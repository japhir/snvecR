## the euler transformation
## :PROPERTIES:
## :CREATED:  [2023-03-24 Fri 15:14]
## :header-args:R: :tangle R/euler.R :comments org :session *R:snvec-R* :exports both :results output :eval no-export
## :END:
## [[file:snvec-3.7.5/snvec-3.7.5.c::=== euler()][euler()]]

#' Euler transformation of spin vector
#'
#' `euler()` transforms a spin vector s from the invariable plane to the
#' instant orbit plane. It uses the formula $s* = A * s$. Setting `inv = TRUE`
#' gives inverse transformation (A^-1 = A' = transpose(A)).
#'
#' @param s The spin vector to be transformed.
#' @param inc  The inclination.
#' @param lan  The long ascending node.
#' @param inv  Invert the output.
#'
#' @returns The transformed vector in the instant orbit plane.
#' @noRd
euler <- function(s, inc, lan, inv = FALSE) {
  ## a <- matrix(
  ##   c(
  ##     cos(lan), sin(lan), 0,
  ##     -cos(inc) * sin(lan), cos(inc) * cos(lan), sin(inc),
  ##     sin(inc) * sin(lan), -sin(inc) * cos(lan), cos(inc)
  ##   ),
  ##   ncol = 3,
  ##   byrow = TRUE
  ## )
  ## if (inv) a <- t(a)
  ## a %*% s
  # this is pretty slow, let's see if we can just use sx sy sx

  a1 <- c(cos(lan), sin(lan), 0)
  a2 <- c(-cos(inc) * sin(lan), cos(inc) * cos(lan), sin(inc))
  a3 <- c(sin(inc) * sin(lan), -sin(inc) * cos(lan), cos(inc))

  # almost direct copy of c-code
  if (inv) {
    # make copy
    b1 <- a1; b2 <- a2; b3 <- a3
    # swap elements
    a1[2] = b2[1]
    a1[3] = b3[1]
    a2[1] = b1[2]
    a2[3] = b3[2]
    a3[1] = b1[3]
    a3[2] = b2[3]
  }

  vp <- c(0.0, 0.0, 0.0)
  for (k in 1:3) {
    vp[1] <- vp[1] + a1[k] * s[k]
    vp[2] <- vp[2] + a2[k] * s[k]
    vp[3] <- vp[3] + a3[k] * s[k]
  }
  vp
  ## c(
  ##   sum(a1 * s[1]),
  ##   sum(a2 * s[2]),
  ##   sum(a3 * s[3])
  ## )

  # TODO: make it take sx, sy, sz separately
  ## c(
  ##   sum(a1 * sx)
  ##   sum(a2 * sy)
  ##   sum(a3 * sz)
  ## )
}
