## for our use-case, the file is
## [[file:snvec-3.7.5/ems-plan3.dat]]

## the top of the file has some lines specifying which columns were used
## 0  7  8  9  12 10 11 15


library(readr)

dat <- read_table("snvec-3.7.5/ems-plan3.dat",
                  comment = "#",
                  skip = 3,
                  col_names = c(
                    "time",              # 0-Time (=Epoch)
                    ## "x1", "x2", "x3", # 1-3
                    ## "v1", "v2", "v3", # 4-6
                    "semimajor_axis",    # 7-SemiMajorAxis
                    "eccentricity",      # 8
                    "inclination",       # 9
                    "long_periapse",     # 12
                    ## "time_periapse",  # 13
                    "long_ascend_node",  # 10
                    "arg_periapse",      # 11
                    ## "peri_distance",  # 12
                    "mean_anomaly"#      # 15
                    ## "true_anomaly",   # 16
                    ## "mean_longitude", # 17
                    ## "true_longitude", # 18
                    ## "mean_latitude",  # 19
                    ## "true_latitude",  # 20
                    ## "mass",           # 21
                    ## "enc_radius",     # 22
                    ## "capt_radius",    # 23
                    ## "id_tag",         # 24
                    ## "jac_index"       # 25
                  ))

head(dat) |> round(2)

## rename some of the names in dat
## :PROPERTIES:
## :CREATED:  [2023-03-24 Fri 14:14]
## :END:
## to make the naming consistent with the C code

dat <- dat |>
  tidylog::rename(
    t  = time,
    aa = semimajor_axis,
    ee = eccentricity,
    inc = inclination,
    lph = long_periapse,
    lan = long_ascend_node,
    arp = arg_periapse,
    mna = mean_anomaly)

## calculate the unwraps for lph and lan
## :PROPERTIES:
## :CREATED:  [2023-03-29 Wed 12:03]
## :END:
## unwrap lph, lan

dat <- dat |>
  dplyr::mutate(lphu = unwrap(lph),
         lanu = unwrap(lan))

## calculate helper parameters
## [[file:snvec-3.7.5/snvec-3.7.5.c::=== fvei()][fvei()]]
## helper parameters as new columns of dat

dat <- dat |>
  dplyr::mutate(age = - t / KY2D, .after = t) |>
  dplyr::mutate(hh = ee * sin(lph / R2D),
         kk = ee * cos(lph / R2D),
         pp = 2 * sin(0.5 * inc / R2D) * sin(lan / R2D),
         qq = 2 * sin(0.5 * inc / R2D) * cos(lan / R2D),
         cc = cos(inc / R2D),
         dd = cos(inc / R2D / 2),
         ## /* nn <- nvec(t): normal to orbit */
         nnx = sin(inc / R2D) * sin(lan / R2D),
         nny = -sin(inc / R2D) * cos(lan / R2D),
         nnz = cos(inc / R2D))

## save the data

ZB18a <- dat
usethis::use_data(ZB18a, overwrite = TRUE, version = 3,
                  ## compress = "bzip2" # 13M
                  ## compress = "gzip" # 14M
                  compress = "xz" # 12M
                  )
