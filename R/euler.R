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
  a <- matrix(
    c(
      cos(lan), sin(lan), 0,
      -cos(inc) * sin(lan), cos(inc) * cos(lan), sin(inc),
      sin(inc) * sin(lan), -sin(inc) * cos(lan), cos(inc)
    ),
    ncol = 3,
    byrow = TRUE
  )
  if (inv) a <- t(a)
  a %*% s
}
