# La11 is available on <http://vo.imcce.fr/insola/earth/online/earth/La2010/La2010c_alkhqp3L.dat>
# the readme says:
# The file La2010X_alkhqp3.dat and La2010X_alkhqp3L.dat contain t, a,l,k,h,q,p
# where t is  the time from J2000 (in kyr)
# a : semi-major axis
# l : mean longitude (expressed in radians)
# k : e  cos (varpi)
# h : e  sin (varpi)
# q : sin(i/2) cos (Omega)
# p : sin(i/2) sin (Omega)

# but note from the README: http://vo.imcce.fr/insola/earth/online/earth/La2010/README.TXT
# The classical elliptical elements a, e, i, l, varpi, Omega
# are taken in the invariant reference frame
# and the origin of time is J2000.

cli::cli_abort(c("i" = "The input OS for snvec must be in the Heliocentric Inertial Reference frame (HCI) (J2000)",
                 "x" = "The La11 solution is in the invariant reference frame"))
# to get eccentricity we'll need a separate file:

# The files La2010X_ecc3.dat and La2010X_ecc3L.dat  contain t, ecc
# where t is  the time from J2000 (in kyr)  and ecc the Earth eccentricity

# The step time is 5kyr in the La2010X_ecc3.dat.
# The step time is 1kyr in the La2010X_ecc3L.dat.


library(readr)
library(snvecR) # for the unwrap function

# this only downloads solution "c".
# Richard told me he prefers b & c the most (<50 Ma?) from Zeebe & Lourens 2019 table 1
dat <- readr::read_table("http://vo.imcce.fr/insola/earth/online/earth/La2010/La2010b_alkhqp3L.dat",
                             col_names = c("time", "a", "l", "k","h", "q", "p"))

ecc <- readr::read_table("http://vo.imcce.fr/insola/earth/online/earth/La2010/La2010b_ecc3L.dat",
                             col_names = c("time", "e"))


dplyr::glimpse(dat)
dplyr::glimpse(ecc)

# back-calculate the input columns
# use wolframalpha to solve this system of equations:
# just e and \varpi:
# https://www.wolframalpha.com/input?i=k%3De%5Ccos%28%5Cvarpi%29%3Bh%3De%5Csin%28%5Cvarpi%29++solve+for+e%2C+%5Cvarpi
# just i and \Omega
# solve \varpi from e and k
## https://www.wolframalpha.com/input?i=k%3De%5Ccos%28%5Cvarpi%29+solve+for+%5Cvarpi+
# solve \varpi from e and h
## https://www.wolframalpha.com/input?i=k%3De%5Ccos%28%5Cvarpi%29+solve+for+%5Cvarpi+
# \varpi = -sin^-1(h/e)+2*pi*n+pi (n is in the set of integers)
# vpi = asin(h/e)+2*pi*n+pi (n is in the set of integers)
# \varpi = sin^-1(h/e)+2*pi*n (n is in the set of integers)
# vpi = asin(h/e)+2*pi*n (n is in the set of integers)
# https://www.wolframalpha.com/input?i=k%3De%5Ccos%28%5Cvarpi%29%3Bh%3De%5Csin%28%5Cvarpi%29%3Bq%3D%5Csin%28i%2F2%29%5Ccos%28%5COmega%29%3Bp%3D%5Csin%28i%2F2%29%5Csin%28%5COmega%29%2C+solve+for+i%2C+%5Cvarpi%2C+%5COmega

n <- 1
La11a <- left_join(dat, ecc) |>
  mutate(
    ## om = from ???
    # wiki says \omega = \arccos(\frac{\vec{n} \dot e}{|n||e|})
    # where n is a vector pointing towards the ascending node and e is eccentricity vector (vector pointing to periapsis)
    om = NA_real_,
    ## oom = from p and q,
    oom = 2 * (pi * n + atan((sqrt(p^2 + q^2) - q) / p)) - 2 * pi,
    ## i = from p and q,
    i = 2 * (2 * pi * n - asin(q / cos(oom)) + pi),
    ## vpi = from e + k and h,
    ## vpi_h = asin(h/e) + 2 * pi * n, #*n (n is in the set of integers)
    ## vpi_k = acos(k/e) + 2 * pi * n, #*n (n is in the set of integers)
    vpi = 2 * pi * n -2 * atan((sqrt(h^2 + k^2) + k) / h) - 2 * pi,
    ## mn = from mean longitude?
    # wiki says: M = E - e sin(E) for which you need the Newton--Rhapson method
    # mean longitude l = \Omega + \omega + M (mean longitude)
    mn = l - oom - om # but we don't have om yet
    # or use true anomaly \nu somehow?
    # \cos\nu = \frac{\cos(E)-e}{1-e\cos(E)}
    # \nu = E + 2*\arctan(\frac{\beta\sin(E)}{1-\beta\cos(E)})
  ) |>
  tidylog::select(time, a, e, i, om, oom, vpi, l, mn, h, k, p, q) |>
  mutate(solution = "La11c")

write_rds(La11a, "~/SurfDrive/Postdoc1/prj/2023-03-23_snvec-R/out/La11.rds")
# TODO: rotate everything so that it aligns with J2000 rather than invariant reference frame

# the stuff below is older!







## rename some of the names in dat
## to make the naming consistent with the C code
dat <- dat |>
  dplyr::left_join(ecc) |>
  tidylog::rename(
    t = t,
    aa = a,
    ee = ecc,
    # not available in this solution
    # which one is the mean longitude?
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
  ## dplyr::mutate(age = -t / KY2D, .after = t) |>
  dplyr::mutate(
    # reverse-engineer lph, lan, inc from khqp
    lph = NA,
    inc = NA,
    lan = NA,
    ## hh = ee * sin(lph / R2D),
    ## kk = ee * cos(lph / R2D),
    ## pp = 2 * sin(0.5 * inc / R2D) * sin(lan / R2D),
    ## qq = 2 * sin(0.5 * inc / R2D) * cos(lan / R2D),
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
