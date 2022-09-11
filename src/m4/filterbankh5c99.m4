# serial 1 uvh5c99.m4
AC_DEFUN([AX_CHECK_FILTERBANKH5C99], [
  AC_PREREQ([2.65])dnl

AC_ARG_WITH([filterbankh5c99],
            AC_HELP_STRING([--with-filterbankh5c99=DIR],
                           [Location of FilterbankH5C99 library]),
            [
              FILTERBANKH5C99DIR="$withval"
              has_filterbankh5=1
            ],
            [
              FILTERBANKH5C99DIR=""
              has_filterbankh5=0
            ])

  if test $has_filterbankh5 = 0; then
    #AC_MSG_RESULT([no])
    AC_MSG_NOTICE([Library FilterbankH5C99 not provided. FilterbankH5C99 will not be linked.])
    filterbankh5_enabled=0;
  else
    # test radiointerferometryc99 before enabling
    AC_CHECK_FILE([${FILTERBANKH5C99DIR}/include/filterbankh5c99.h],
                  # Found
                  AC_SUBST(FILTERBANKH5C99_INCDIR,${FILTERBANKH5C99DIR}/include),
                  AC_MSG_ERROR([filterbankh5c99.h header file not found])
    )
                  

    orig_LDFLAGS="${LDFLAGS}"
    LDFLAGS="${orig_LDFLAGS} -L${FILTERBANKH5C99DIR}/"
    AC_CHECK_LIB([filterbankh5c99], [filterbankh5_open],
                # Found
                AC_SUBST(FILTERBANKH5C99_LIBDIR,${FILTERBANKH5C99DIR}/),
                # Not found there, check under lib
                AS_UNSET(ac_cv_lib_filterbankh5c99_filterbankh5_open)
                LDFLAGS="${orig_LDFLAGS} -L${FILTERBANKH5C99DIR}/lib"
                AC_CHECK_LIB([filterbankh5c99], [filterbankh5_open],
                  # Found
                  AC_SUBST(FILTERBANKH5C99_LIBDIR,${FILTERBANKH5C99DIR}/lib),
                  # Not found there, try under lib/x86_64-linux-gnu
                  AS_UNSET(ac_cv_lib_filterbankh5c99_filterbankh5_open)
                  LDFLAGS="${orig_LDFLAGS} -L${FILTERBANKH5C99DIR}/lib/x86_64-linux-gnu"
                  AC_CHECK_LIB([filterbankh5c99], [filterbankh5_open],
                    # Found
                    AC_SUBST(FILTERBANKH5C99_LIBDIR,${FILTERBANKH5C99DIR}/lib/x86_64-linux-gnu),
                    # Not found there
                    AS_UNSET(ac_cv_lib_filterbankh5c99_filterbankh5_open)
                    AC_MSG_ERROR([FILTERBANKH5C99_LIBDIR library not found])
                  )
                )
    )
    LDFLAGS="${orig_LDFLAGS}"
    filterbankh5_enabled=1;
  fi

  AS_IF([test $filterbankh5_enabled = 1],
  [
    AM_CONDITIONAL(FILTERBANKH5C99_ENABLED, true)
    AC_DEFINE(FILTERBANKH5C99_ENABLED,[],[Use FilterbankH5C99])
  ],
  [
    AM_CONDITIONAL(FILTERBANKH5C99_ENABLED, false)
  ])
])
