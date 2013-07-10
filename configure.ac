AC_INIT(GenAlgPLSConfig, version-1.2)

AC_PROG_CXX
AC_LANG(C++)

AC_CHECK_HEADERS(stdint.h inttypes.h limits.h climits)

oldCFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -D__STDC_LIMIT_MACROS"

# Check if unsigned long long is supported (its not standard yet)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  unsigned long long x = 100LL; return 0; }]])
	],
	[
		AC_MSG_RESULT([Type unsigned long long is supported])
		AC_DEFINE_UNQUOTED(HAVE_UNSIGNED_LONG_LONG)
	],
	[AC_MSG_RESULT([Type unsigned long long is not supported])]
)

# Check if support for uint16_t and uint8_t is complete
if test "$ac_cv_header_stdint_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_PROGRAM([[#include <stdint.h>]],
					[[uint16_t x = UINT16_MAX; uint8_t y = UINT8_MAX;]])
	],
	[
		AC_MSG_RESULT([uint16_t and unit8_t are available in stdint.h])
		AC_DEFINE_UNQUOTED(HAVE_UINT8_16_MAX)
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint16_t and uint8_t are not available in stdint.h])]
	)
fi
if test "$ac_cv_header_inttypes_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_PROGRAM([[#include <inttypes.h>]],
					[[uint16_t x = UINT16_MAX; uint8_t y = UINT8_MAX;]])
	],
	[
		AC_MSG_RESULT([uint16_t and unit8_t are available in inttypes.h])
		AC_DEFINE_UNQUOTED(HAVE_UINT8_16_MAX)
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint16_t and uint8_t are not available in inttypes.h])]
	)
fi

if test "$suitable_integer" != "yes"; then
	AC_MSG_FAILURE([no suitable integral type available])
fi;

CFLAGS="$oldCFLAGS"

# Check for builtin function to "count trailing zeros" (ctz)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_ctzl(2UL); return 0; }]])
	],
	[
		AC_MSG_RESULT([Found builtin CTZ method (for unsigned long)])
		AC_DEFINE_UNQUOTED(HAVE_GCC_CTZL)
	],
	[AC_MSG_RESULT([No builtin CTZ method (for unsigned long)])]
)

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_ctzll(2ULL); return 0; }]])
	],
	[
		AC_MSG_RESULT([Found builtin CTZ method (for unsigned long long)])
		AC_DEFINE_UNQUOTED(HAVE_GCC_CTZLL)
	],
	[AC_MSG_RESULT([No builtin CTZ method (for unsigned long long)])]
)

AC_CHECK_FUNCS(ffsl)
AC_CHECK_FUNCS(ffsll)

# Check for builtin function to "count set bits" (popcountl)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_popcountl(10UL); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin popcount method (for unsigned long)])
		AC_DEFINE_UNQUOTED(HAVE_BUILTIN_POPCOUNTL)
	],
	[AC_MSG_RESULT([No builtin popcount method (for unsigned long)])]
)

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_popcountll(10ULL); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin popcount method (for unsigned long long)])
		AC_DEFINE_UNQUOTED(HAVE_BUILTIN_POPCOUNTLL)
	],
	[AC_MSG_RESULT([No builtin popcount method (for unsigned long long)])]
)

cp confdefs.h src/autoconfig.h