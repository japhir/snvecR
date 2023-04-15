
<!-- README.md is generated from README.Rmd. Please edit that file -->

# snvecR

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/snvecR)](https://CRAN.R-project.org/package=snvecR)
[![GPL-3](https://img.shields.io/github/license/japhir/snvecR?logo=gnu&.svg)](https://github.com/japhir/snvecR/blob/master/LICENSE.md)
[![release](https://img.shields.io/github/v/release/japhir/snvecR.svg)](https://github.com/japhir/snvecR/releases)
[![R-CMD-check](https://github.com/japhir/snvecR/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/japhir/snvecR/actions/workflows/check-standard.yaml)
[![Launch
binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/japhir/snvecR/main)
<!-- badges: end -->

Easily calculate precession and obliquity from an orbital solution (OS,
defaults to ZB18a from Zeebe and Lourens (2019)) and assumed or
reconstructed values for tidal dissipation (T<sub>d</sub>) and dynamical
ellipticity (E<sub>d</sub>). This is a translation and adaptation of the
C-code in the supplementary material to Zeebe and Lourens (2022), with
furher details on the methodology described in Zeebe (2022).

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
#> • `tend` = -1000 ka
#> • `ed` = 1
#> • `td` = 0
#> • `orbital_solution` = "ZB18a"
#> • `tres` = 0.4 kyr
#> • `tolerance` = 1e-07
#> ℹ started at "2023-04-14 16:26:22"
#> Final values:
#> • s[1][2][3]: 0.404197400723194 -0.0537088738295803 0.91303387030935
#> • s-error = |s|-1: -5.44863786333671e-05
#> Final values:
#> • obliquity: 0.413056573207875 rad
#> • precession: -0.562236553023642 rad
#> ℹ stopped at "2023-04-14 16:26:25"
#> ℹ total duration: 2.76s
```

see `?snvec` for further documentation.

Here we create a quick plot of the calculated climatic precession with
the eccentricity envelope:

``` r
library(ggplot2)
solution |>
  ggplot(aes(x = age, y = cp)) +
  labs(x = "Age (ka)", y = "(-)", colour = "Orbital Element") +
  # the age scale goes from old to young
  scale_x_reverse() +
  # plot climatic precession
  geom_line(aes(colour = "Climatic Precession")) +
  # add the (interpolated) eccentricity envelope
  geom_line(aes(y = eei, colour = "Eccentricity")) +
  scale_color_discrete(type = c("skyblue", "black")) +
  theme(legend.pos = c(.9, .95))
```

<img src="man/figures/README-plot-1.png" width="100%" />

# References

Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
Paleocene–Eocene boundary age constrained by geology and astronomy.
*Science*, 365(6456), 926–929. <https://doi.org/10.1126/science.aax0612>

Zeebe, R. E., & Lourens, L. J. (2022). A deep-time dating tool for
paleo-applications utilizing obliquity and precession cycles: The role
of dynamical ellipticity and tidal dissipation. *Paleoceanography and
Paleoclimatology*, e2021PA004349. <https://doi.org/10.1029/2021PA004349>

Zeebe, R. E. (2022). Reduced Variations in Earth’s and Mars’ Orbital
Inclination and Earth’s Obliquity from 58 to 48 Myr ago due to Solar
System Chaos. *The Astronomical Journal*, 164(3), 107.
<https://doi.org/10.3847/1538-3881/ac80f8>

Wikipedia page on Orbital Elements:
<https://en.wikipedia.org/wiki/Orbital_elements>
