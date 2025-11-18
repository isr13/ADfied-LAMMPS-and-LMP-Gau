
  const char features[] = {"\n"
"CXX_FEATURE:"
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 500 && __cplusplus >= 201402L
"1"
#else
"0"
#endif
"cxx_aggregate_default_initializers\n"
"CXX_FEATURE:"
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 407 && __cplusplus >= 201103L
"1"
#else
"0"
#endif
"cxx_alias_templates\n"
"CXX_FEATURE:"
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L
"1"
#else
"0"
#endif
"cxx_alignas\n"
"CXX_FEATURE:"
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L
"1"
#else
"0"
#endif
"cxx_alignof\n"
"CXX_FEATURE:"
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L
"1"
#else
"0"
#endif
"cxx_attributes\n"
"CXX_FEATURE:"
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 409 && __cplusplus > 201103L
"1"
#else
"0"
#endif
"cxx_attribute_deprecated\n"
"CXX_FEATURE:"
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 404 && (__cplusplus >= 201103L || (defined(__GXX_EXPERIMENTAL_CXX0X__) && __GXX_EXPERIMENTAL_CXX0X__))
"1"
#else
"0"
#endif
"cxx_auto_type\n"
