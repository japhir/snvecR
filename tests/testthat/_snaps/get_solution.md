# get_solution() works

    Code
      head(get_solution(orbital_solution = "ZB18a"))
    Output
      # A tibble: 6 x 20
              t   age    aa     ee   inc   lph   lan   arp   mna  lphu  lanu      hh
          <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
      1       0   0    1.00 0.0167  7.15  27.3  180. -153. -2.45  27.3  180. 0.00767
      2 -146100   0.4  1.00 0.0169  7.15  26.1 -180. -154.  1.27  26.1  180. 0.00742
      3 -292200   0.8  1.00 0.0171  7.14  24.7 -179. -156.  5.22  24.7  181. 0.00713
      4 -438300   1.2  1.00 0.0172  7.12  23.7 -179. -158.  8.75  23.7  181. 0.00690
      5 -584400   1.6  1.00 0.0173  7.11  22.1 -178. -160. 12.8   22.1  182. 0.00653
      6 -730500   2    1.00 0.0175  7.10  21.0 -178. -161. 16.4   21.0  182. 0.00627
      # i 8 more variables: kk <dbl>, pp <dbl>, qq <dbl>, cc <dbl>, dd <dbl>,
      #   nnx <dbl>, nny <dbl>, nnz <dbl>

