LIBTOOL=libtool
RANLIB=ranlib

noinst_HEADERS = eqlm550a.h filters.h rms.h db.h blo.h pitchscale.h waveguide_nl.h biquad.h iir.h buffer.h ls_filter.h

noinst_LIBRARIES = libeqlm550a.a libfilters.a librms.a libdb.a libblo.a libpitchscale.a libiir.a

librms_a_SOURCES = rms.c
libdb_a_SOURCES = db.c
libblo_a_SOURCES = blo.c
libfilters_a_SOURCES = filters.c
libeqlm550a_a_SOURCES = eqlm550a.c
libblo_a_CFLAGS = -I@top_srcdir@
libpitchscale_a_SOURCES = pitchscale.c
libpitchscale_a_CFLAGS = $(FFTW_CFLAGS)
#libpitchscale_a_LIBS = $(FFTW_LIBS)
libiir_a_SOURCES = iir.c
libiir_a_CPPFLAGS = -I@top_srcdir@
#libbuffer_a_SOURCES = buffer.c
# Disable autoheader. 
AUTOHEADER=echo
