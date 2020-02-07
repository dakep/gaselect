This release updates the autoconf script to remove the deprecated call to `R CMD config CXXCPP`,
addresses an issue with potentially uninitialized memory (reported by valgrind),
and updates the `cleanup` script to use `/bin/sh` instead of `/bin/bash`.

## Test environments

* local OS X 10.14.6, R 3.6.1
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

