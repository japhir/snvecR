# get_solution() works

    Code
      head(get_solution(astronomical_solution = "PT-ZB18a", quiet = TRUE))
    Output
      # A tibble: 6 x 20
              t t_kyr    aa     ee   inc   lph   lan   arp   mna  lphu  lanu      hh     kk          pp     qq    cc    dd         nnx   nny   nnz
          <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>  <dbl>       <dbl>  <dbl> <dbl> <dbl>       <dbl> <dbl> <dbl>
      1       0   0    1.00 0.0167  7.15  27.3  180. -153. -2.45  27.3  180. 0.00767 0.0148  0.00000164 -0.125 0.992 0.998  0.00000163 0.125 0.992
      2 -146100  -0.4  1.00 0.0169  7.15  26.1 -180. -154.  1.27  26.1  180. 0.00742 0.0151 -0.000902   -0.125 0.992 0.998 -0.000900   0.124 0.992
      3 -292200  -0.8  1.00 0.0171  7.14  24.7 -179. -156.  5.22  24.7  181. 0.00713 0.0155 -0.00180    -0.124 0.992 0.998 -0.00180    0.124 0.992
      4 -438300  -1.2  1.00 0.0172  7.12  23.7 -179. -158.  8.75  23.7  181. 0.00690 0.0157 -0.00270    -0.124 0.992 0.998 -0.00270    0.124 0.992
      5 -584400  -1.6  1.00 0.0173  7.11  22.1 -178. -160. 12.8   22.1  182. 0.00653 0.0161 -0.00359    -0.124 0.992 0.998 -0.00359    0.124 0.992
      6 -730500  -2    1.00 0.0175  7.10  21.0 -178. -161. 16.4   21.0  182. 0.00627 0.0163 -0.00449    -0.124 0.992 0.998 -0.00448    0.123 0.992

---

    Code
      head(get_solution(astronomical_solution = "PT-ZB18a", quiet = FALSE, force = TRUE))
    Message
      i The astronomical solution PT-ZB18a has not been downloaded.
      i Reading 'PT-ZB18a.dat' from website <http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat>.
      i Calculating helper columns.
      i The cache directory is 'transformed-for-CI'.
      i Saved 'PT-ZB18a.dat' to cache.
      i Saved cleaned-up 'PT-ZB18a.csv' to cache.
      > Saved solution with helper columns 'PT-ZB18a.rds' to cache.
    Output
      # A tibble: 6 x 20
              t t_kyr    aa     ee   inc   lph   lan   arp   mna  lphu  lanu      hh     kk          pp     qq    cc    dd         nnx   nny   nnz
          <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>  <dbl>       <dbl>  <dbl> <dbl> <dbl>       <dbl> <dbl> <dbl>
      1       0   0    1.00 0.0167  7.15  27.3  180. -153. -2.45  27.3  180. 0.00767 0.0148  0.00000164 -0.125 0.992 0.998  0.00000163 0.125 0.992
      2 -146100  -0.4  1.00 0.0169  7.15  26.1 -180. -154.  1.27  26.1  180. 0.00742 0.0151 -0.000902   -0.125 0.992 0.998 -0.000900   0.124 0.992
      3 -292200  -0.8  1.00 0.0171  7.14  24.7 -179. -156.  5.22  24.7  181. 0.00713 0.0155 -0.00180    -0.124 0.992 0.998 -0.00180    0.124 0.992
      4 -438300  -1.2  1.00 0.0172  7.12  23.7 -179. -158.  8.75  23.7  181. 0.00690 0.0157 -0.00270    -0.124 0.992 0.998 -0.00270    0.124 0.992
      5 -584400  -1.6  1.00 0.0173  7.11  22.1 -178. -160. 12.8   22.1  182. 0.00653 0.0161 -0.00359    -0.124 0.992 0.998 -0.00359    0.124 0.992
      6 -730500  -2    1.00 0.0175  7.10  21.0 -178. -161. 16.4   21.0  182. 0.00627 0.0163 -0.00449    -0.124 0.992 0.998 -0.00448    0.123 0.992

