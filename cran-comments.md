This release fixes the issue with the incorrect C++ flags in the configure script as reported by Prof. Ripley.
The package now complies with C++17.

## Test environments

* macOS 12.6.2, R 4.2.2 Patched (2023-02-03 r83757)
* win-builder (devel and release)
* Rhub
  * Debian Linux, R-release, GCC
  * Debian Linux, R-devel, GCC
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE.

* on some systems the size of the shared library is > 5MB.

