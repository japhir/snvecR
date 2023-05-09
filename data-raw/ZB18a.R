## for our use-case, the file is
## [[file:snvec-3.7.5/ems-plan3.dat]]
# also available on <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat>

## the top of the file has some lines specifying which columns were used
## 0  7  8  9  12 10 11 15
# where the numbers correspond to HNBody column indicators

library(readr)
library(snvecR) # for the unwrap function

dat <- read_table("http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat",
                  # also locally at
                  #"snvec-3.7.5/ems-plan3.dat"
  comment = "#",
  skip = 3, # skip those column indicators
  # set the column names manually to match the header comments
  col_names = c(
    "Time",              # 0-Time (=Epoch)
    ## "x1", "x2", "x3", # 1-3
    ## "v1", "v2", "v3", # 4-6
    "SemiMajorAxis",     # 7
    "Eccentricity",      # 8
    "Inclination",       # 9
    "LongPeriapse",      # 12
    ## "TimePeriapse",   # 13
    "LongAscendNode",    # 10
    "ArgPeriapse",       # 11
    ## "PeriDistance",   # 12
    "MeanAnomaly"        # 15
    ## "TrueAnomaly",    # 16
    ## "MeanLongitude",  # 17
    ## "TrueLongitude",  # 18
    ## "MeanLatitude",   # 19
    ## "TrueLatitude",   # 20
    ## "Mass"            # 21
    ## "EncRadius",      # 22
    ## "CaptRadius",     # 23
    ## "IdTag",          # 24
    ## "JacIndex"        # 25
  )
)

dplyr::glimpse(dat)

## rename some of the names in dat
## to make the naming consistent with the C code
dat <- dat |>
  tidylog::rename(
    t = Time,
    aa = SemiMajorAxis,
    ee = Eccentricity,
    inc = Inclination,
    lph = LongPeriapse,
    lan = LongAscendNode,
    arp = ArgPeriapse,
    mna = MeanAnomaly
  )

dplyr::glimpse(dat)

## calculate the unwraps for lph and lan
dat <- dat |>
  dplyr::mutate(
    # this function is not exported so we need to use three colons
    lphu = snvecR:::unwrap(lph),
    lanu = snvecR:::unwrap(lan)
  )

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
ZB18a <- dat

usethis::use_data(ZB18a,
  overwrite = TRUE, version = 3,
  ## compress = "bzip2" # 13M
  ## compress = "gzip" # 14M
  compress = "xz" # 12M
)

rm(dat)

# all of this is now a part of the fuction get_ZB18a()!
