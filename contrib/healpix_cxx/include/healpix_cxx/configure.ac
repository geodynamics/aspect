AC_INIT([healpix_cxx],[3.82.0])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AM_MAINTAINER_MODE([enable])

dnl
dnl Needed for linking on Windows.
dnl Protect with m4_ifdef because AM_PROG_AR is required in
dnl autoconf >= 1.12 when using -Wall, but the macro is
dnl absent in old versions of autoconf.
dnl
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

LT_INIT
AC_CONFIG_MACRO_DIR([m4])

dnl
dnl By default, install the headers into a subdirectory of
dnl ${prefix}/include to avoid possible header filename collisions.
dnl
includedir="${includedir}/${PACKAGE_NAME}"

dnl
dnl Enable silent build rules if this version of Automake supports them
dnl
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

m4_version_prereq(2.70, [AC_PROG_CC], [AC_PROG_CC_C99])
AC_OPENMP
CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
CFLAGS="$CFLAGS $OPENMP_CFLAGS"

AC_PROG_CXX
AC_LANG_PUSH([C++])
AX_CXX_COMPILE_STDCXX_11
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
AC_LANG_POP([C++])

PKG_CHECK_MODULES([ZLIB], [zlib], , [ZLIB_LIBS=-lz])
CFLAGS="$CFLAGS $ZLIB_CFLAGS"
LIBS="$LIBS $ZLIB_LIBS"

dnl
dnl Determine the compiler flags (CFITSIO_CFLAGS)
dnl and linker flags (CFITSIO_LIBS) required to build against libcfitsio.
dnl
dnl Add help messages to describe the variables CFITSIO_CFLAGS and CFITSIO_LIBS,
dnl which can be provided as environment variables or as command-line arguments.
dnl These are output variables of PKG_CHECK_MODULES, but they can be overridden.
dnl
AC_ARG_VAR([CFITSIO_CFLAGS], [compiler flags needed for building with cfitsio, if not found by pkg-config])
AC_ARG_VAR([CFITSIO_LIBS], [linker flags needed for building with cfitsio, if not found by pkg-config])
dnl
dnl Now check for cfitsio using pkg-config. If cfitsio is not found with
dnl pkg-config, CFITSIO_LIBS defaults to "-lcfitsio".
dnl
PKG_CHECK_MODULES([CFITSIO], [cfitsio], , [CFITSIO_LIBS=-lcfitsio])
dnl
dnl Linking against libcfitsio requires libm too.
dnl
CFITSIO_LIBS="$CFITSIO_LIBS -lm"
dnl
dnl With pkg-config < 0.24, PKG_CHECK_MODULES does not AC_SUBST the _CFLAGS
dnl and _LIBS variables.
dnl
AC_SUBST([CFITSIO_CFLAGS])
AC_SUBST([CFITSIO_LIBS])
dnl
dnl Perform a test compile and link to determine if cfitsio is usable. If
dnl pkg-config did not find cfitsio, then the user may still have provided
dnl CFITSIO_CFLAGS and CFITSIO_LIBS manually. If pkg-config did find cfitsio,
dnl it does not hurt to test that cfitsio is usable.
dnl
save_CFLAGS=$CFLAGS
save_LIBS=$LIBS
CFLAGS="$CFLAGS $CFITSIO_CFLAGS"
LIBS="$LIBS $CFITSIO_LIBS"
AC_CHECK_LIB([cfitsio], [ffgnrwll], [:], AC_MSG_ERROR([could not find the cfitsio library]))
AC_CHECK_HEADERS([fitsio.h], [:], AC_MSG_ERROR([could not find the cfitsio header file]))
CFLAGS=$save_CFLAGS
LIBS=$save_LIBS

dnl
dnl Determine the compiler flags (SHARP_CFLAGS)
dnl and linker flags (SHARP_LIBS) required to build against libcfitsio.
dnl
dnl Add help messages to describe the variables SHARP_CFLAGS and SHARP_LIBS,
dnl which can be provided as environment variables or as command-line arguments.
dnl These are output variables of PKG_CHECK_MODULES, but they can be overridden.
dnl
AC_ARG_VAR([SHARP_CFLAGS], [compiler flags needed for building with libsharp, if not found by pkg-config])
AC_ARG_VAR([SHARP_LIBS], [linker flags needed for building with libsharp, if not found by pkg-config])
dnl
dnl Now check for libsharp using pkg-config. If libsharp is not found with
dnl pkg-config, SHARP_LIBS defaults to "-lsharp".
dnl
PKG_CHECK_MODULES([SHARP], [libsharp], , [SHARP_LIBS=-lsharp])
dnl
dnl Linking against libsharp requires libm too.
dnl
SHARP_LIBS="$SHARP_LIBS -lm"
dnl
dnl With pkg-config < 0.24, PKG_CHECK_MODULES does not AC_SUBST the _CFLAGS
dnl and _LIBS variables.
dnl
AC_SUBST([SHARP_CFLAGS])
AC_SUBST([SHARP_LIBS])
dnl
dnl Perform a test compile and link to determine if libsharp is usable. If
dnl pkg-config did not find libsharp, then the user may still have provided
dnl SHARP_CFLAGS and SHARP_LIBS manually. If pkg-config did find libsharp,
dnl it does not hurt to test that libsharp is usable.
dnl
save_CFLAGS=$CFLAGS
save_LIBS=$LIBS
CFLAGS="$CFLAGS $SHARP_CFLAGS"
LIBS="$LIBS $SHARP_LIBS"
AC_CHECK_LIB([sharp], [sharp_execute], [:], AC_MSG_ERROR([could not find the libsharp library]))
AC_CHECK_HEADERS([libsharp/sharp.h], [:], AC_MSG_ERROR([could not find the libsharp header file]))
CFLAGS=$save_CFLAGS
LIBS=$save_LIBS

LT_INIT

dnl
dnl Create pkgconfig .pc file.
dnl
AX_CREATE_PKGCONFIG_INFO(,,,,[$OPENMP_CFLAGS $SHARP_CFLAGS $CFITSIO_CFLAGS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
