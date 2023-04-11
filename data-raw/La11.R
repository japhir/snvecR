# La11 is available on <http://vo.imcce.fr/insola/earth/online/earth/La2010/La2010d_alkhqp3L.dat>
# the readme says:
# The file La2010X_alkhqp3.dat and La2010X_alkhqp3L.dat contain t, a,l,k,h,q,p
# where t is  the time from J2000 (in kyr)
# a : semi-major axis
# l : mean longitude (expressed in radians)
# k : e  cos (varpi)
# h : e  sin (varpi)
# q : sin(i/2) cos (Omega)
# p : sin(i/2) sin (Omega)

library(readr)
library(snvecR) # for the unwrap function

dat <- readr::read_table("http://vo.imcce.fr/insola/earth/online/earth/La2010/La2010d_alkhqp3L.dat",
                             col_names = c("t", "a", "l", "k","h", "q", "p"))

dplyr::glimpse(dat)

## rename some of the names in dat
## to make the naming consistent with the C code
dat <- dat |>
  tidylog::rename(
    t = t,
    aa = a,
    ## ee = e, # NOT AVAILABLE in this one!
    ## inc = Inclination,
    ## lph = LongPeriapse,
    ## lan = LongAscendNode,
    ## arp = ArgPeriapse,
    ## mna = MeanAnomaly
    kk = k, hh = h, qq = q, pp = p)

dplyr::glimpse(dat)

## calculate the unwraps for lph and lan
## dat <- dat |>
##   dplyr::mutate(
##     # this function is not exported so we need to use three colons
##     lphu = snvecR:::unwrap(lph),
##     lanu = snvecR:::unwrap(lan)
##   )

## calculate helper parameters
## [[file:snvec-3.7.5/snvec-3.7.5.c::=== fvei()][fvei()]]
## helper parameters as new columns of dat
dat <- dat |>
  dplyr::mutate(age = -t / KY2D, .after = t) |>
  dplyr::mutate(
    hh = ee * sin(lph / R2D),
    kk = ee * cos(lph / R2D),
    pp = 2 * sin(0.5 * inc / R2D) * sin(lan / R2D),
    qq = 2 * sin(0.5 * inc / R2D) * cos(lan / R2D),
    cc = cos(inc / R2D),
    dd = cos(inc / R2D / 2),
    ## /* nn <- nvec(t): normal to orbit */
    nnx = sin(inc / R2D) * sin(lan / R2D),
    nny = -sin(inc / R2D) * cos(lan / R2D),
    nnz = cos(inc / R2D)
  )

## save the data
# we want the object to have the correct name
La11 <- dat

usethis::use_data(La11,
  overwrite = TRUE, version = 3,
  ## compress = "bzip2" # 13M
  ## compress = "gzip" # 14M
  compress = "xz" # 12M
)

rm(dat)
