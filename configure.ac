AC_INIT(GenAlgPLSConfig, version-1.0)

AC_PROG_CXX
AC_LANG(C++)

AC_CHECK_HEADER([stdint.h])
AC_CHECK_HEADER([inttypes.h])

oldCFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -D__STDC_LIMIT_MACROS"

# Check if support for uint64_t is complete
if test "$ac_cv_header_stdint_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_SOURCE([[#include <stdint.h> int main() {  uint64_t x = UINT64_MAX; return 0; }]])
	],
	[
		AC_MSG_RESULT([uint64_t is available in inttypes.h])
		AC_DEFINE_UNQUOTED(HAVE_UINT64_MAX)
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint64_t is not available in stdint.h])],
	)
fi
if test "$ac_cv_header_inttypes_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_SOURCE([[#include <inttypes.h> int main() {  uint64_t x = UINT64_MAX; return 0; }]])
	],
	[
		AC_MSG_RESULT([uint64_t is available in inttypes.h])
		AC_DEFINE_UNQUOTED(HAVE_UINT64_MAX)
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint64_t is not available in inttypes.h])],
	)
fi

# Check if support for uint64_t is complete
if test "$ac_cv_header_stdint_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_SOURCE([[#include <stdint.h> int main() {  uint32_t x = UINT32_MAX; return 0; }]])
	],
	[
		AC_MSG_RESULT([uint32_t is available in stdint.h])
		AC_DEFINE_UNQUOTED(HAVE_UINT32_MAX)
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint32_t is not available in stdint.h])],
	)
fi
if test "$ac_cv_header_inttypes_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_SOURCE([[#include <inttypes.h> int main() {  uint32_t x = UINT32_MAX; return 0; }]])
	],
	[
		AC_MSG_RESULT([uint32_t is available in stdint.h])
		AC_DEFINE_UNQUOTED(HAVE_UINT32_MAX)
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint32_t is not available in inttypes.h])],
	)
fi

if test "$suitable_integer" != "yes"; then
	AC_MSG_FAILURE([no suitable integral type available])
fi;

CFLAGS="$oldCFLAGS"

# Check for builtin function to "count trailing zeros" (ctz)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_ctz(2), y = __builtin_ctzl(2); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin CTZ method])
		AC_DEFINE_UNQUOTED(HAVE_GCC_CTZ)
	],
	[AC_MSG_RESULT([No builtin CTZ method])])

AC_CHECK_FUNCS(ffs)
AC_CHECK_FUNCS(ffsl)


# Check for builtin popcount

# Check for builtin function to "count trailing zeros" (ctz)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_popcountl(10); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin popcount method])
		AC_DEFINE_UNQUOTED(HAVE_BUILTIN_POPCOUNT)
	],
	[AC_MSG_RESULT([No builtin popcount method])])

cp confdefs.h src/autoconfig.h