# snvecR basic call works

    Code
      print(dplyr::select(snvec(tend = -49, ed = 1, td = 0,
        tres = 1, quiet = TRUE, output = "all"), dplyr::all_of(
        c("age", "sx", "sy", "sz", "epl", "phi", "cp"))), n = 50)
    Output
      # A tibble: 50 x 7
           age      sx         sy    sz   epl     phi        cp
         <dbl>   <dbl>      <dbl> <dbl> <dbl>   <dbl>     <dbl>
       1     0  0.385   0.212     0.898 0.409  0       0.0163  
       2     1  0.350   0.302     0.887 0.411  0.243   0.0171  
       3     2  0.293   0.381     0.877 0.414  0.486   0.0163  
       4     3  0.219   0.445     0.868 0.416  0.727   0.0140  
       5     4  0.132   0.489     0.862 0.418  0.967   0.0105  
       6     5  0.0360  0.511     0.859 0.419  1.21    0.00573 
       7     6 -0.0629  0.511     0.857 0.421  1.45    0.000463
       8     7 -0.159   0.487     0.859 0.422  1.68   -0.00499 
       9     8 -0.247   0.441     0.863 0.423  1.92   -0.0101  
      10     9 -0.322   0.376     0.869 0.423  2.16   -0.0145  
      11    10 -0.379   0.295     0.877 0.423  2.40   -0.0176  
      12    11 -0.416   0.203     0.887 0.422  2.63   -0.0193  
      13    12 -0.429   0.105     0.897 0.422  2.87   -0.0194  
      14    13 -0.420   0.00696   0.908 0.420  3.11   -0.0180  
      15    14 -0.387  -0.0858    0.918 0.419 -2.94   -0.0150  
      16    15 -0.334  -0.168     0.928 0.417 -2.70   -0.0108  
      17    16 -0.262  -0.235     0.936 0.415 -2.46   -0.00569 
      18    17 -0.177  -0.283     0.943 0.412 -2.22   -0.000231
      19    18 -0.0839 -0.309     0.947 0.409 -1.97    0.00523 
      20    19  0.0128 -0.312     0.950 0.407 -1.73    0.0102  
      21    20  0.107  -0.291     0.951 0.404 -1.49    0.0142  
      22    21  0.193  -0.249     0.949 0.401 -1.24    0.0170  
      23    22  0.265  -0.187     0.946 0.398 -0.993   0.0184  
      24    23  0.320  -0.110     0.941 0.395 -0.744   0.0182  
      25    24  0.354  -0.0228    0.935 0.393 -0.494   0.0166  
      26    25  0.365   0.0701    0.928 0.391 -0.242   0.0136  
      27    26  0.352   0.162     0.922 0.390  0.0104  0.00966 
      28    27  0.317   0.248     0.915 0.388  0.264   0.00517 
      29    28  0.261   0.323     0.910 0.388  0.517   0.000383
      30    29  0.189   0.381     0.905 0.388  0.771  -0.00423 
      31    30  0.104   0.419     0.902 0.388  1.03   -0.00830 
      32    31  0.0122  0.434     0.901 0.389  1.28   -0.0116  
      33    32 -0.0807  0.426     0.901 0.390  1.53   -0.0139  
      34    33 -0.169   0.395     0.903 0.392  1.78   -0.0151  
      35    34 -0.247   0.343     0.906 0.394  2.03   -0.0152  
      36    35 -0.310   0.273     0.911 0.397  2.28   -0.0143  
      37    36 -0.354   0.189     0.916 0.400  2.53   -0.0125  
      38    37 -0.376   0.0959    0.922 0.403  2.78   -0.0100  
      39    38 -0.375   0.0000166 0.927 0.406  3.03   -0.00704 
      40    39 -0.350  -0.0933    0.932 0.409 -3.01   -0.00375 
      41    40 -0.304  -0.178     0.936 0.412 -2.77   -0.000399
      42    41 -0.237  -0.250     0.939 0.415 -2.53    0.00287 
      43    42 -0.156  -0.305     0.940 0.418 -2.29    0.00595 
      44    43 -0.0629 -0.338     0.939 0.420 -2.05    0.00861 
      45    44  0.0355 -0.349     0.936 0.422 -1.81    0.0108  
      46    45  0.134  -0.337     0.932 0.424 -1.58    0.0125  
      47    46  0.227  -0.302     0.926 0.425 -1.34    0.0136  
      48    47  0.309  -0.245     0.919 0.426 -1.10    0.0140  
      49    48  0.376  -0.171     0.911 0.426 -0.868   0.0137  
      50    49  0.424  -0.0839    0.902 0.426 -0.632   0.0126  
