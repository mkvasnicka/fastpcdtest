# fastpcdtest

A package for fast and robust calculation of Pesaran's CD test in R.
It substitutes `pcdtest` in **plm** package in two situations:

1. when the panel is huges (since it is faster than pcdtest in plm package) and
2. when the panel is highly unbalanced (since is calculates CD test even when
    some observations do not have overlap in their time dimension).

The package provides only one function, `pcdtest`.

Beware: this is just an alpha version. It may contain bugs and
the interface may change.
