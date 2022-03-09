As requested by Prof. Brian Ripley this release updates the C++ preprocessor flags for RcppArmadillo to fix the warning
```
WARNING: option ARMA_DONT_USE_CXX11 ignored [-W#pragma-messages]
```

CRAN checks on `r-devel-linux-x86_64-debian-gcc` raise several warnings:

* Warnings related to a mismatch between allocator and delete have been addressed by replacing the corresponding raw pointers with std::vector.
* Warnings regarding the use of deprecated C++11 features are due to Rcpp. Rcpp developers are aware of these warnings and are working on a fix.

## Test environments

* macOS 12.2.1, R 4.1.2
* win-builder (devel and release)
* Rhub
  * Debian Linux, R-release, GCC
  * Debian Linux, R-devel, GCC
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE.

* on some systems the size of the shared library is > 5MB.

