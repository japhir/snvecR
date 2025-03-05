# get_solution() can return a dataframe

    Code
      get_solution(astronomical_solution = ZB18a_head, quiet = TRUE)
    Output
      # A tibble: 2 x 20
              t  time    aa     ee   inc   lph   lan   arp   mna  lphu  lanu      hh
          <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
      1       0   0    1.00 0.0167  7.15  27.3  180. -153. -2.45  27.3  180. 0.00767
      2 -146100  -0.4  1.00 0.0169  7.15  26.1 -180. -154.  1.27  26.1  180. 0.00742
      # i 8 more variables: kk <dbl>, pp <dbl>, qq <dbl>, cc <dbl>, dd <dbl>,
      #   nnx <dbl>, nny <dbl>, nnz <dbl>

# get_solution() can load eccentricity solutions

    Code
      head(get_solution(astronomical_solution = "ZB17a", quiet = TRUE, force = TRUE))
    Output
      # A tibble: 6 x 3
         time    ecc   inc
        <dbl>  <dbl> <dbl>
      1   0   0.0167  7.15
      2  -1.6 0.0173  7.11
      3  -3.2 0.0179  7.04
      4  -4.8 0.0184  6.95
      5  -6.4 0.0188  6.84
      6  -8   0.0192  6.70

---

    Code
      head(get_solution(astronomical_solution = "ZB18a-100", quiet = TRUE, force = TRUE))
    Output
      # A tibble: 6 x 3
         time    ecc   inc
        <dbl>  <dbl> <dbl>
      1   0   0.0167  7.15
      2  -1.6 0.0173  7.11
      3  -3.2 0.0179  7.04
      4  -4.8 0.0184  6.95
      5  -6.4 0.0188  6.84
      6  -8   0.0192  6.70

---

    Code
      head(get_solution(astronomical_solution = "ZB18a-300", quiet = TRUE, force = TRUE))
    Output
      # A tibble: 6 x 3
         time    ecc   inc
        <dbl>  <dbl> <dbl>
      1   0   0.0167  7.15
      2  -1.6 0.0173  7.11
      3  -3.2 0.0179  7.04
      4  -4.8 0.0184  6.95
      5  -6.4 0.0188  6.84
      6  -8   0.0192  6.70

---

    Code
      head(get_solution(astronomical_solution = "ZB20a", quiet = TRUE, force = TRUE))
    Output
      # A tibble: 6 x 3
         time    ecc   inc
        <dbl>  <dbl> <dbl>
      1   0   0.0167  7.15
      2  -1.6 0.0173  7.11
      3  -3.2 0.0179  7.04
      4  -4.8 0.0184  6.95
      5  -6.4 0.0188  6.84
      6  -8   0.0192  6.70

---

    Code
      head(get_solution(astronomical_solution = "La10d", quiet = TRUE, force = TRUE))
    Message
      i Relying on astrochron to get solution "La10d"
      i We do not cache these results.
      ! astrochron converts time from -kyr to ka by default.
      i Output has column names "Time_ka" and "ecc_LA10d"
    Output
      # A tibble: 6 x 2
        Time_ka ecc_LA10d
          <dbl>     <dbl>
      1       0    0.0167
      2       1    0.0172
      3       2    0.0175
      4       3    0.0178
      5       4    0.0182
      6       5    0.0185

# get_solution() can load full solutions

    Code
      head(get_solution(astronomical_solution = "full-ZB18a", quiet = FALSE, force = TRUE))
    Message
      i The astronomical solution "full-ZB18a" has not been cached.
      i Reading 'full-ZB18a.dat' from website <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat>.
      i Calculating helper columns.
      i The cache directory is 'transformed-for-CI'.
      i Saved astronomical solution with helper columns 'full-ZB18a.rds' to cache.
      i Future calls to `get_solution("full-ZB18a")` will read from the cache.
      ! If you want to read from scratch, specify `force = TRUE`.
    Output
      # A tibble: 6 x 20
              t  time    aa     ee   inc   lph   lan   arp   mna  lphu  lanu      hh     kk          pp     qq    cc    dd         nnx   nny   nnz
          <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>  <dbl>       <dbl>  <dbl> <dbl> <dbl>       <dbl> <dbl> <dbl>
      1       0   0    1.00 0.0167  7.15  27.3  180. -153. -2.45  27.3  180. 0.00767 0.0148  0.00000164 -0.125 0.992 0.998  0.00000163 0.125 0.992
      2 -146100  -0.4  1.00 0.0169  7.15  26.1 -180. -154.  1.27  26.1  180. 0.00742 0.0151 -0.000902   -0.125 0.992 0.998 -0.000900   0.124 0.992
      3 -292200  -0.8  1.00 0.0171  7.14  24.7 -179. -156.  5.22  24.7  181. 0.00713 0.0155 -0.00180    -0.124 0.992 0.998 -0.00180    0.124 0.992
      4 -438300  -1.2  1.00 0.0172  7.12  23.7 -179. -158.  8.75  23.7  181. 0.00690 0.0157 -0.00270    -0.124 0.992 0.998 -0.00270    0.124 0.992
      5 -584400  -1.6  1.00 0.0173  7.11  22.1 -178. -160. 12.8   22.1  182. 0.00653 0.0161 -0.00359    -0.124 0.992 0.998 -0.00359    0.124 0.992
      6 -730500  -2    1.00 0.0175  7.10  21.0 -178. -161. 16.4   21.0  182. 0.00627 0.0163 -0.00449    -0.124 0.992 0.998 -0.00448    0.123 0.992

---

    Code
      head(get_solution(astronomical_solution = "full-ZB18a", quiet = TRUE))
    Output
      # A tibble: 6 x 20
              t  time    aa     ee   inc   lph   lan   arp   mna  lphu  lanu      hh     kk          pp     qq    cc    dd         nnx   nny   nnz
          <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>  <dbl>       <dbl>  <dbl> <dbl> <dbl>       <dbl> <dbl> <dbl>
      1       0   0    1.00 0.0167  7.15  27.3  180. -153. -2.45  27.3  180. 0.00767 0.0148  0.00000164 -0.125 0.992 0.998  0.00000163 0.125 0.992
      2 -146100  -0.4  1.00 0.0169  7.15  26.1 -180. -154.  1.27  26.1  180. 0.00742 0.0151 -0.000902   -0.125 0.992 0.998 -0.000900   0.124 0.992
      3 -292200  -0.8  1.00 0.0171  7.14  24.7 -179. -156.  5.22  24.7  181. 0.00713 0.0155 -0.00180    -0.124 0.992 0.998 -0.00180    0.124 0.992
      4 -438300  -1.2  1.00 0.0172  7.12  23.7 -179. -158.  8.75  23.7  181. 0.00690 0.0157 -0.00270    -0.124 0.992 0.998 -0.00270    0.124 0.992
      5 -584400  -1.6  1.00 0.0173  7.11  22.1 -178. -160. 12.8   22.1  182. 0.00653 0.0161 -0.00359    -0.124 0.992 0.998 -0.00359    0.124 0.992
      6 -730500  -2    1.00 0.0175  7.10  21.0 -178. -161. 16.4   21.0  182. 0.00627 0.0163 -0.00449    -0.124 0.992 0.998 -0.00448    0.123 0.992

# get_solution() can load PT solutions

    Code
      head(get_solution(astronomical_solution = "PT-ZB18a(1,1)", quiet = FALSE, force = TRUE))
    Message
      i The astronomical solution "PT-ZB18a(1.0000,1.0000)" has not been cached.
      i Reading 'PT-ZB18a(1.0000,1.0000).dat' from website <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/ZB18a/asc/PT.De1.0000Td1.0000.dat>.
      i The cache directory is 'transformed-for-CI'.
      i Saved astronomical solution with helper columns 'PT-ZB18a(1.0000,1.0000).rds' to cache.
      i Future calls to `get_solution("PT-ZB18a(1.0000,1.0000)")` will read from the cache.
      ! If you want to read from scratch, specify `force = TRUE`.
    Output
      # A tibble: 6 x 4
          time   epl    phi     cp
         <dbl> <dbl>  <dbl>  <dbl>
      1  0     0.409 0      0.0163
      2 -0.379 0.410 0.0924 0.0167
      3 -0.760 0.411 0.185  0.0170
      4 -1.14  0.412 0.277  0.0171
      5 -1.52  0.413 0.369  0.0169
      6 -1.90  0.413 0.461  0.0165

---

    Code
      head(get_solution(astronomical_solution = "PT-ZB18a(1,1)", quiet = TRUE))
    Output
      # A tibble: 6 x 4
          time   epl    phi     cp
         <dbl> <dbl>  <dbl>  <dbl>
      1  0     0.409 0      0.0163
      2 -0.379 0.410 0.0924 0.0167
      3 -0.760 0.411 0.185  0.0170
      4 -1.14  0.412 0.277  0.0171
      5 -1.52  0.413 0.369  0.0169
      6 -1.90  0.413 0.461  0.0165

# get_solution() can load ZB23.Rxx solutions

    Code
      head(get_solution(astronomical_solution = "ZB23.R01", quiet = FALSE, force = TRUE))
    Message
      i The astronomical solution "ZB23.R01" has not been cached.
      i Reading 'ZB23.R01.dat' from website <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/3.5Gyr/ZB23-N64-eiop/ZB23.R01.eiop.dat.zip>.
      i Downloading any of the ZB23.RXX solutions will take some time.
      i Zip files are about 154 MB.
    Output
      Continue downloading and caching? (Yes/no/cancel) 
    Message
      i The cache directory is 'transformed-for-CI'.
      i Saved astronomical solution with helper columns 'ZB23.R01.rds' to cache.
      i Future calls to `get_solution("ZB23.R01")` will read from the cache.
      ! If you want to read from scratch, specify `force = TRUE`.
    Output
      # A tibble: 6 x 5
         time    ecc   inc   epl     cp
        <dbl>  <dbl> <dbl> <dbl>  <dbl>
      1   0   0.0167 0.125 0.409 0.0163
      2  -0.4 0.0169 0.125 0.410 0.0168
      3  -0.8 0.0170 0.124 0.411 0.0170
      4  -1.2 0.0172 0.124 0.412 0.0170
      5  -1.6 0.0173 0.124 0.413 0.0168
      6  -2   0.0175 0.124 0.414 0.0163

---

    Code
      head(get_solution(astronomical_solution = "ZB23.R01", quiet = TRUE))
    Output
      # A tibble: 6 x 5
         time    ecc   inc   epl     cp
        <dbl>  <dbl> <dbl> <dbl>  <dbl>
      1   0   0.0167 0.125 0.409 0.0163
      2  -0.4 0.0169 0.125 0.410 0.0168
      3  -0.8 0.0170 0.124 0.411 0.0170
      4  -1.2 0.0172 0.124 0.412 0.0170
      5  -1.6 0.0173 0.124 0.413 0.0168
      6  -2   0.0175 0.124 0.414 0.0163

