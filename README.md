
<!-- README.md is generated from README.Rmd. Please edit that file -->

# snvecR

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/snvecR)](https://CRAN.R-project.org/package=snvecR)
<!-- badges: end -->

The goal of snvecR is to calculate precession and obliquity from an
orbital solution (defaults to ZB18a) and assumed or reconstructed values
for tidal dissipation (td) and dynamical ellipticity (ed).

## Installation

You can install the development version of snvecR like so:

``` r
remotes::install_github("japhir/snvecR")
```

## Example

Here’s the main function that does the work in action:

``` r
library(snvecR)
solution <- snvec()
#> Integration parameters:
#> • `tend` = -1000 kyr
#> • `ed` = 1
#> • `td` = 0
#> • `orbital_solution` = "ZB18a"
#> • `tres` = 0.4 kyr
#> • `tolerance` = 1e-07
#> ℹ started at "2023-04-12 09:05:04"
#> Final values:
#> • s[1][2][3]: 0.404197400723194 -0.0537088738295803 0.91303387030935
#> • s-error = |s|-1: -5.44863786333671e-05
#> Final values:
#> • obliquity: 0.413056573207875 rad
#> • precession: -0.562236553023642 rad
#> ℹ stopped at "2023-04-12 09:05:05"
#> ℹ total duration: 1.19s
```

see `?snvec` for further documentation.

Here we create a quick plot of the calculated climatic precession with
the eccentricity envelope:

``` r
library(ggplot2)
solution |>
  ggplot(aes(x = age, y = cp)) +
  # the age scale goes from old to young
  scale_x_reverse() +
  # plot climatic precession
  geom_line() +
  # add the (interpolated) eccentricity envelope
  geom_line(aes(y = eei))
```

<img src="man/figures/README-plot-1.png" width="100%" />
