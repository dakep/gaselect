This release fixes the warnings generated on r-devel-linux-x86_64-debian-gcc and r-devel-linux-x86_64-debian-clang 
related to deprecated C++ features.
Thanks to Kurt Hornik for alerting me to these warnings.

## Test environments

* macOS 12.3, R 4.1.3 Patched (2022-03-10 r82100)
* win-builder (devel and release)
* Rhub
  * Debian Linux, R-release, GCC
  * Debian Linux, R-devel, GCC
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE.

* on some systems the size of the shared library is > 5MB.

