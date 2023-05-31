# snvecR 3.7.6
* Remove cache directory after last example (even though it's in a donttest environment).

# snvecR 3.7.5
* Added a `NEWS.md` file to track changes to the package.
* Note: version number is consistent with snvec C package.
* Implemented first public version based on C-code in [Zeebe and Lourens
  2022](https://doi.org/10.1029/2021PA004349), as explained .
* Added `snvec()` function that does all the work.
* Added `ZB18a` dataset with orbital solution.
* Added tests using `testthat`.
* Added `snvec()` `output` parameter with choice to select `"ode"`, `"all"`, or `"nice"` (default).
* Added a vignette with a grid of variations on Td and Ed.
* Added a bookdown website rendering.
* Released a version to Zenodo and assigned a doi.
* Removed `ZB18a` dataset from the package because it was too large.
* Added functions get_solution, get_ZB18a, prepare_solution.
* Added caching code for the orbital solution.
* Remove the cache directory after running tests, so that reproducible environments remain clean.
