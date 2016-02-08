#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT(GSATools, 1.0.0, [alessandro.pandini@gmail.com])
AC_CONFIG_SRCDIR([src/g_sa_encode.c])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CC_C_O
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m], [main], [], 
             AC_MSG_ERROR([Could not find required C math library.]))

AC_ARG_WITH(gsl, [  --with-gsl=<path>      prefix of GSL installation],
                 [  CPPFLAGS="$CPPFLAGS -I $withval/include"
                    LDFLAGS="$LDFLAGS -L $withval/lib"])

AC_CHECK_HEADER([gsl/gsl_math.h], [],
             AC_MSG_ERROR([Could not find required GSL header files.
                           Consider using --with-gsl=<my-path-to-gsl>
                           or installing the GSL if required.
                           Read the INSTALL file for details.]))

AC_CHECK_LIB([gslcblas], [main], [], 
            AC_MSG_ERROR([Could not find required GSL BLAS library.
                          Consider using --with-gsl=<my-path-to-gsl>
                          or installing the GSL if required.
                          Read the INSTALL file for details.]))

AC_CHECK_LIB([gsl], [main], [],
             AC_MSG_ERROR([Could not find required GSL library.
                          Consider using --with-gsl=<my-path-to-gsl>
                          or installing the GSL if required.
                          Read the INSTALL file for details.]))

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([strcasecmp])

# check GROMACS installation dir.
AC_ARG_VAR(GMXDIR, [path to gromacs installation])
AC_MSG_CHECKING([whether GMXDIR has been set])
if test -z "$GMXDIR" ; then
    GMXDIR="$GMXDIR"
fi
if test -d "$GMXDIR" ; then
	AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR( [  GMXDIR is undefined but required for compilation.
                     Set GMXDIR according to your Gromacs path.
                     Read the INSTALL file for details.] )
fi

# set the default compilation to O3
CFLAGS="-O3 -fomit-frame-pointer -finline-functions"

# enable debug
AC_ARG_ENABLE(debug, [  --enable-debug          Turn on debugging],
                [ac_use_debug="yes"], [ac_use_debug="no"])
if test "$ac_use_debug" = "yes"; then
   CFLAGS="-ggdb"
   AC_DEFINE(DEBUG,1,[Compiling with DEBUG settings.])
fi

# enable profiling
AC_ARG_ENABLE(profiling, [  --enable-profiling      Turn on profiling],
                [ac_use_profiling="yes"], [ac_use_profiling="no"])
if test "$ac_use_profiling" = "yes"; then
   CFLAGS="-pg -O2"
   AC_DEFINE(PROFILING,1,[Compiling with PROFILING settings.])
fi

# MPI commpiler search with ax_mpi.m4 
# http://www.gnu.org/software/autoconf-archive/ax_mpi.html
sinclude(ax_mpi.m4)

# enable MPI
AC_ARG_ENABLE(mpi, [  --enable-mpi         Turn on MPI],
                [ac_use_mpi="yes"], [ac_use_mpi="no"])
AC_ARG_WITH(mpi, [  --with-mpi=<path>      prefix of MPI installation],
                 [  MPICC="$withval/bin/mpicc"
                    CPPFLAGS="$CPPFLAGS -I $withval/include"
                    LDFLAGS="$LDFLAGS -L $withval/lib"])

if test "$ac_use_mpi" = "yes"; then
   AX_MPI([], AC_MSG_ERROR([Could not find required MPI library.
                         Consider using --with-mpi=<my-path-to-mpi>.]))
   CC="$MPICC"
   LIBS="$MPILIBS $LIBS"
   AC_DEFINE(MPI,1,[Compiling with MPI settings.])
   AC_DEFINE(GMX_LIB_MPI,1,[Using system installed version of MPI.])
fi

# Output files
AC_CONFIG_FILES(Makefile src/Makefile tests/Makefile scripts/Makefile)
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

