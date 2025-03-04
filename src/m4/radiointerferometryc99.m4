# serial 1 uvh5c99.m4
AC_DEFUN([AX_CHECK_RADIOINTERFEROMETRY99], [
  AC_PREREQ([2.65])dnl

AC_ARG_WITH([radiointerferometry],
            AC_HELP_STRING([--with-radiointerferometry=DIR],
                           [Location of RadioInterferometryC99 library]),
            [
              RADIOINTERFEROMETRYDIR="$withval"
              has_radiointerferometry=1
            ],
            [
              RADIOINTERFEROMETRYDIR=""
              has_radiointerferometry=0
            ])

  radiointerferometry_enabled=0;  

  if test $has_radiointerferometry = 0; then
    if test $uvh5c99_enabled = 1; then
      
      # test radiointerferometryc99 before enabling
      AC_CHECK_FILE([${UVH5C99_INCDIR}/radiointerferometryc99.h],
                    # Found
                    AC_SUBST(RADIOINTERFEROMETRYC99_INCDIR,${UVH5C99_INCDIR})
                    orig_LDFLAGS="${LDFLAGS}"
                    LDFLAGS="${orig_LDFLAGS} -L${UVH5C99_LIBDIR}"
                    AC_CHECK_LIB([radiointerferometryc99], [calc_independent_astrom],
                                # Found
                                AC_SUBST(RADIOINTERFEROMETRYC99_LIBDIR,${UVH5C99_LIBDIR})
                                radiointerferometry_enabled=1,
                                # Not found there
                                AC_MSG_NOTICE([RADIOINTERFEROMETRYC99_LIBDIR library not found alongside UVH5 library])
                                radiointerferometry_enabled=0
                    )
                    LDFLAGS="${orig_LDFLAGS}"                    
                    ,
                    AC_MSG_NOTICE([RADIOINTERFEROMETRYC99_INCDIR header not found alongside UVH5 header])
                    radiointerferometry_enabled=0
      )

      #radiointerferometryc99 is a subproject of uvh5c99
      RADIOINTERFEROMETRYDIR="${UVH5C99_LIBDIR}/subprojects/radiointerferometryc99"
      has_radiointerferometry=1;
    fi
  fi

  if test $has_radiointerferometry = 0; then
    #AC_MSG_RESULT([no])
    AC_MSG_NOTICE([Library RadioInterferometryC99 not provided. RadioInterferometryC99 will not be linked.])
  else
    if test $radiointerferometry_enabled = 0; then
      # test radiointerferometryc99 before enabling
      AC_CHECK_FILE([${RADIOINTERFEROMETRYDIR}/include/radiointerferometryc99.h],
                    # Found
                    AC_SUBST(RADIOINTERFEROMETRYC99_INCDIR,${RADIOINTERFEROMETRYDIR}/include),
                    # Not found there, check UVH5C99_LIBDIR/../subprojects/radiointerferometryc99/include
                    AC_CHECK_FILE([${UVH5C99_LIBDIR}/../subprojects/radiointerferometryc99/include/radiointerferometryc99.h],
                                  # Found
                                  AC_SUBST(RADIOINTERFEROMETRYC99_INCDIR,${UVH5C99_LIBDIR}/../subprojects/radiointerferometryc99/include),
                                  # Not found there, error
                                  AC_MSG_ERROR([radiointerferometryc99.h header file not found])))

      orig_LDFLAGS="${LDFLAGS}"
      LDFLAGS="${orig_LDFLAGS} -L${RADIOINTERFEROMETRYDIR}/lib/x86_64-linux-gnu/"
      AC_CHECK_LIB([radiointerferometryc99], [calc_independent_astrom],
                  # Found
                  AC_SUBST(RADIOINTERFEROMETRYC99_LIBDIR,${RADIOINTERFEROMETRYDIR}/lib/x86_64-linux-gnu/),
                  # Not found there, check RADIOINTERFEROMETRYDIR
                  AS_UNSET(ac_cv_lib_radiointerferometryc99_calc_independent_astrom)
                  LDFLAGS="${orig_LDFLAGS} -L${RADIOINTERFEROMETRYDIR}"
                  AC_CHECK_LIB([radiointerferometryc99], [calc_independent_astrom],
                              # Found
                              AC_SUBST(RADIOINTERFEROMETRYC99_LIBDIR,${RADIOINTERFEROMETRYDIR}),
                              # Not found there, error
                              AC_MSG_ERROR([RADIOINTERFEROMETRYC99_LIBDIR library not found])))
      LDFLAGS="${orig_LDFLAGS}"
      radiointerferometry_enabled=1;
    fi
  fi

  AS_IF([test $radiointerferometry_enabled = 1],
  [
    AM_CONDITIONAL(RADIOINTERFEROMETRYC99_ENABLED, true)
    AC_DEFINE(RADIOINTERFEROMETRYC99_ENABLED,[],[Use RadioInterferometryC99])
  ],
  [
    AM_CONDITIONAL(RADIOINTERFEROMETRYC99_ENABLED, false)
  ])
])
