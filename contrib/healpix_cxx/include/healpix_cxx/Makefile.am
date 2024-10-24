ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libhealpix_cxx.la

src_cxxsupport = \
  cxxsupport/announce.cc \
  cxxsupport/geom_utils.cc \
  cxxsupport/string_utils.cc \
  cxxsupport/ls_image.cc \
  cxxsupport/paramfile.cc \
  cxxsupport/pointing.cc \
  cxxsupport/rotmatrix.cc \
  cxxsupport/trafos.cc \
  cxxsupport/walltimer.cc \
  cxxsupport/wigner.cc \
  cxxsupport/error_handling.cc \
  cxxsupport/alloc_utils.h \
  cxxsupport/announce.h \
  cxxsupport/arr.h \
  cxxsupport/bstream.h \
  cxxsupport/datatypes.h \
  cxxsupport/error_handling.h \
  cxxsupport/geom_utils.h \
  cxxsupport/levels_facilities.h \
  cxxsupport/ls_image.h \
  cxxsupport/lsconstants.h \
  cxxsupport/math_utils.h \
  cxxsupport/openmp_support.h \
  cxxsupport/paramfile.h \
  cxxsupport/planck_rng.h \
  cxxsupport/pointing.h \
  cxxsupport/rangeset.h \
  cxxsupport/crangeset.h \
  cxxsupport/rotmatrix.h \
  cxxsupport/safe_cast.h \
  cxxsupport/share_utils.h \
  cxxsupport/sse_utils_cxx.h \
  cxxsupport/string_utils.h \
  cxxsupport/trafos.h \
  cxxsupport/vec3.h \
  cxxsupport/walltimer.h \
  cxxsupport/wigner.h \
  cxxsupport/xcomplex.h \
  cxxsupport/compress_utils.h \
  cxxsupport/colour.h \
  cxxsupport/linear_map.h \
  cxxsupport/sort_utils.h
src_cxxsupport_fits = \
  cxxsupport/fitshandle.cc \
  cxxsupport/fitshandle.h
src_healpix_cxx= \
  Healpix_cxx/alm.cc \
  Healpix_cxx/alm.h \
  Healpix_cxx/alm_healpix_tools.cc \
  Healpix_cxx/alm_healpix_tools.h \
  Healpix_cxx/alm_powspec_tools.cc \
  Healpix_cxx/alm_powspec_tools.h \
  Healpix_cxx/healpix_tables.cc \
  Healpix_cxx/healpix_tables.h \
  Healpix_cxx/healpix_base.cc \
  Healpix_cxx/healpix_base.h \
  Healpix_cxx/healpix_map.cc \
  Healpix_cxx/healpix_map.h \
  Healpix_cxx/powspec.cc \
  Healpix_cxx/powspec.h \
  Healpix_cxx/moc_query.h \
  Healpix_cxx/moc_query.cc \
  Healpix_cxx/weight_utils.h \
  Healpix_cxx/weight_utils.cc \
  Healpix_cxx/mask_tools.cc
src_healpix_cxx_fits= \
  Healpix_cxx/healpix_data_io.cc \
  Healpix_cxx/healpix_data_io.h \
  Healpix_cxx/healpix_map_fitsio.cc \
  Healpix_cxx/healpix_map_fitsio.h \
  Healpix_cxx/alm_fitsio.cc \
  Healpix_cxx/alm_fitsio.h \
  Healpix_cxx/powspec_fitsio.cc \
  Healpix_cxx/powspec_fitsio.h \
  Healpix_cxx/moc_fitsio.h \
  Healpix_cxx/moc_fitsio.cc

EXTRA_DIST = \
  cxxsupport/font_data.inc

include_HEADERS = \
  cxxsupport/arr.h \
  cxxsupport/xcomplex.h \
  cxxsupport/datatypes.h \
  cxxsupport/pointing.h \
  cxxsupport/geom_utils.h \
  cxxsupport/rangeset.h \
  cxxsupport/fitshandle.h \
  cxxsupport/alloc_utils.h \
  cxxsupport/compress_utils.h \
  cxxsupport/math_utils.h \
  cxxsupport/error_handling.h \
  cxxsupport/lsconstants.h \
  cxxsupport/safe_cast.h \
  cxxsupport/vec3.h \
  cxxsupport/wigner.h \
  cxxsupport/rotmatrix.h \
  cxxsupport/sse_utils_cxx.h \
  cxxsupport/string_utils.h \
  cxxsupport/planck_rng.h \
  Healpix_cxx/alm.h \
  Healpix_cxx/alm_fitsio.h \
  Healpix_cxx/alm_healpix_tools.h \
  Healpix_cxx/alm_powspec_tools.h \
  Healpix_cxx/healpix_base.h \
  Healpix_cxx/healpix_data_io.h \
  Healpix_cxx/healpix_map.h \
  Healpix_cxx/healpix_map_fitsio.h \
  Healpix_cxx/healpix_tables.h \
  Healpix_cxx/powspec.h \
  Healpix_cxx/powspec_fitsio.h \
  Healpix_cxx/moc.h \
  Healpix_cxx/moc_fitsio.h \
  Healpix_cxx/moc_query.h \
  Healpix_cxx/weight_utils.h \
  Healpix_cxx/mask_tools.h \
  Healpix_cxx/alice3.h

libhealpix_cxx_la_SOURCES = $(src_cxxsupport) $(src_cxxsupport_fits) $(src_healpix_cxx) $(src_healpix_cxx_fits) $(src_healpix_cxx_mod)

libhealpix_cxx_la_LIBADD = $(SHARP_LIBS) $(CFITSIO_LIBS)

# format is "current:revision:age"
# - Always increase the revision value.
# - Increase the current value whenever an interface has been added, removed or changed.
# - Increase the age value only if the changes made to the ABI are backward compatible.
libhealpix_cxx_la_LDFLAGS = -version-info 3:4:0

bin_PROGRAMS = alm2map_cxx anafast_cxx calc_powspec hotspots_cxx map2tga \
  median_filter_cxx mult_alm rotalm_cxx smoothing_cxx syn_alm_cxx udgrade_cxx \
  udgrade_harmonic_cxx compute_weights needlet_tool alice3
check_PROGRAMS = hpxtest
alm2map_cxx_SOURCES = Healpix_cxx/alm2map_cxx.cc Healpix_cxx/alm2map_cxx_module.cc
alm2map_cxx_LDADD = libhealpix_cxx.la
anafast_cxx_SOURCES = Healpix_cxx/anafast_cxx.cc Healpix_cxx/anafast_cxx_module.cc
anafast_cxx_LDADD = libhealpix_cxx.la
calc_powspec_SOURCES = Healpix_cxx/calc_powspec.cc Healpix_cxx/calc_powspec_module.cc
calc_powspec_LDADD = libhealpix_cxx.la
hotspots_cxx_SOURCES = Healpix_cxx/hotspots_cxx.cc Healpix_cxx/hotspots_cxx_module.cc
hotspots_cxx_LDADD = libhealpix_cxx.la
map2tga_SOURCES = Healpix_cxx/map2tga.cc Healpix_cxx/map2tga_module.cc
map2tga_LDADD = libhealpix_cxx.la
median_filter_cxx_SOURCES = Healpix_cxx/median_filter_cxx.cc Healpix_cxx/median_filter_cxx_module.cc
median_filter_cxx_LDADD = libhealpix_cxx.la
mult_alm_SOURCES = Healpix_cxx/mult_alm.cc Healpix_cxx/mult_alm_module.cc
mult_alm_LDADD = libhealpix_cxx.la
rotalm_cxx_SOURCES = Healpix_cxx/rotalm_cxx.cc
rotalm_cxx_LDADD = libhealpix_cxx.la
smoothing_cxx_SOURCES = Healpix_cxx/smoothing_cxx.cc Healpix_cxx/smoothing_cxx_module.cc
smoothing_cxx_LDADD = libhealpix_cxx.la
syn_alm_cxx_SOURCES = Healpix_cxx/syn_alm_cxx.cc Healpix_cxx/syn_alm_cxx_module.cc
syn_alm_cxx_LDADD = libhealpix_cxx.la
udgrade_cxx_SOURCES = Healpix_cxx/udgrade_cxx.cc Healpix_cxx/udgrade_cxx_module.cc
udgrade_cxx_LDADD = libhealpix_cxx.la
udgrade_harmonic_cxx_SOURCES = Healpix_cxx/udgrade_harmonic_cxx.cc Healpix_cxx/udgrade_harmonic_cxx_module.cc
udgrade_harmonic_cxx_LDADD = libhealpix_cxx.la
compute_weights_SOURCES = Healpix_cxx/compute_weights.cc Healpix_cxx/compute_weights_module.cc
compute_weights_LDADD = libhealpix_cxx.la
needlet_tool_SOURCES = Healpix_cxx/needlet_tool.cc Healpix_cxx/needlet_tool_module.cc
needlet_tool_LDADD = libhealpix_cxx.la
alice3_SOURCES = Healpix_cxx/alice3.cc
alice3_LDADD = libhealpix_cxx.la
hpxtest_SOURCES=Healpix_cxx/hpxtest.cc
hpxtest_LDADD = libhealpix_cxx.la

TESTS=hpxtest

AM_CPPFLAGS = $(CFITSIO_CFLAGS) $(SHARP_CFLAGS) -I$(top_srcdir) -I$(top_srcdir)/cxxsupport

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

doc:
	cd $(top_srcdir)/docsrc && for i in cxxsupport Healpix_cxx; do doxygen $${i}.dox; mv htmldoc $${i}; done

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
