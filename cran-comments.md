As requested by Prof. Brian Ripley this release updates the C++ preprocessor flags for RcppArmadillo to fix the warning
```
WARNING: option ARMA_DONT_USE_CXX11 ignored [-W#pragma-messages]
```

## Test environments

* macOS 12.2.1, R 4.1.2
* win-builder (devel and release)
* Rhub
  * Debian Linux, R-release, GCC
  * Debian Linux, R-devel, GCC ASAN/UBSAN
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-patched, 32/64 bit
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE.

* on some systems the size of the shared library is > 5MB.

