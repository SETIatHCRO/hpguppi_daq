# serial 1 filterbankc99.m4
AC_DEFUN([AX_CHECK_FILTERBANKC99], [
  AC_PREREQ([2.65])dnl

AC_ARG_WITH([filterbankc99],
            AC_HELP_STRING([--with-filterbankc99=DIR],
                           [Location of FilterbankC99 library]),
            [
              FILTERBANKC99DIR="$withval"
              has_filterbank=1
            ],
            [
              FILTERBANKC99DIR=""
              has_filterbank=0
            ])

  if test $has_filterbank = 0; then
    #AC_MSG_RESULT([no])
    AC_MSG_NOTICE([Library FilterbankC99 not provided. FilterbankC99 will not be linked.])
    filterbank_enabled=0;
  else
    # test filterbankc99 before enabling
    AC_CHECK_FILE([${FILTERBANKC99DIR}/include/filterbankc99.h],
                  # Found
                  AC_SUBST(FILTERBANKC99_INCDIR,${FILTERBANKC99DIR}/include),
                  AC_MSG_ERROR([filterbankc99.h header file not found])
    )
                  

    orig_LDFLAGS="${LDFLAGS}"
    LDFLAGS="${orig_LDFLAGS} -lh5dsc99 -L${FILTERBANKC99DIR}/"
    AC_CHECK_LIB([filterbankc99], [filterbank_h5_open_explicit],
                # Found
                AC_SUBST(FILTERBANKC99_LIBDIR,${FILTERBANKC99DIR}/),
                # Not found there, check under lib
                AS_UNSET(ac_cv_lib_filterbankc99_filterbank_h5_open_explicit)
                LDFLAGS="${orig_LDFLAGS} -lh5dsc99 -L${FILTERBANKC99DIR}/lib"
                AC_CHECK_LIB([filterbankc99], [filterbank_h5_open_explicit],
                  # Found
                  AC_SUBST(FILTERBANKC99_LIBDIR,${FILTERBANKC99DIR}/lib),
                  # Not found there, try under lib/x86_64-linux-gnu
                  AS_UNSET(ac_cv_lib_filterbankc99_filterbank_h5_open_explicit)
                  LDFLAGS="${orig_LDFLAGS} -lh5dsc99 -L${FILTERBANKC99DIR}/lib/x86_64-linux-gnu"
                  AC_CHECK_LIB([filterbankc99], [filterbank_h5_open_explicit],
                    # Found
                    AC_SUBST(FILTERBANKC99_LIBDIR,${FILTERBANKC99DIR}/lib/x86_64-linux-gnu),
                    # Not found there
                    AS_UNSET(ac_cv_lib_filterbankc99_filterbank_h5_open_explicit)
                    AC_MSG_ERROR([FILTERBANKC99_LIBDIR library not found])
                  )
                )
    )
    LDFLAGS="${orig_LDFLAGS}"
    filterbank_enabled=1;
  fi

  AS_IF([test $filterbank_enabled = 1],
  [
    AM_CONDITIONAL(FILTERBANKC99_ENABLED, true)
    AC_DEFINE(FILTERBANKC99_ENABLED,[],[Use FilterbankC99])
  ],
  [
    AM_CONDITIONAL(FILTERBANKC99_ENABLED, false)
  ])
])
