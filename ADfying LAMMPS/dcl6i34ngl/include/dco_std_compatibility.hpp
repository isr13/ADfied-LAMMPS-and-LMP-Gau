/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#pragma once

#ifdef DCO_STD_COMPATIBILITY
namespace std {
    using dco::internal::ceil;
    using dco::internal::floor;

    using dco::internal::min;
    using dco::internal::max;

    using dco::internal::sin;
    using dco::internal::cos;
    using dco::internal::tan;
    using dco::internal::cosh;
    using dco::internal::sinh;
    using dco::internal::asin;
    using dco::internal::acos;
    using dco::internal::exp;
    using dco::internal::atan;
    using dco::internal::tanh;
    using dco::internal::sqrt;
    using dco::internal::log;
#ifndef _WIN32
    using dco::internal::isfinite;
    using dco::internal::isnan;
    using dco::internal::isnormal;
    using dco::internal::isinf;

    using dco::internal::erf;
    using dco::internal::erfc;
    using dco::internal::expm1;
    using dco::internal::asinh;
    using dco::internal::acosh;
    using dco::internal::atanh;
    using dco::internal::log1p;
    using dco::internal::log10;
#endif
    using dco::internal::fabs;
    using dco::internal::abs;
    using dco::internal::atan2;
    using dco::internal::pow;
    using dco::internal::hypot;

#ifdef DCO_CPP11
#define DCOi_CONSTEXPR constexpr
#else
#define DCOi_CONSTEXPR
#endif

    template<class DCO_TAPE_REAL, class DATA_HANDLER>
    struct numeric_limits<dco::internal::active_type<DCO_TAPE_REAL, DATA_HANDLER> > : public numeric_limits<DCO_TAPE_REAL> {
        typedef dco::internal::active_type<DCO_TAPE_REAL, DATA_HANDLER> DCO_T;
        static DCOi_CONSTEXPR DCO_T min() { return numeric_limits<typename DCO_T::VALUE_TYPE>::min(); }
        static DCOi_CONSTEXPR DCO_T max() { return numeric_limits<typename DCO_T::VALUE_TYPE>::max(); }
        static DCOi_CONSTEXPR DCO_T epsilon() { return numeric_limits<typename DCO_T::VALUE_TYPE>::epsilon() ; }
        static DCOi_CONSTEXPR DCO_T round_error() { return numeric_limits<typename DCO_T::VALUE_TYPE>::round_error() ; }
        static DCOi_CONSTEXPR DCO_T infinity() { return numeric_limits<typename DCO_T::VALUE_TYPE>::infinity(); }
        static DCOi_CONSTEXPR DCO_T quiet_NaN() { return numeric_limits<typename DCO_T::VALUE_TYPE>::quiet_NaN(); }
        static DCOi_CONSTEXPR DCO_T signaling_NaN() { return numeric_limits<typename DCO_T::VALUE_TYPE>::signaling_NaN(); }
        static DCOi_CONSTEXPR DCO_T denorm_min() { return numeric_limits<typename DCO_T::VALUE_TYPE>::denorm_min(); }
#ifdef DCO_CPP11
        static DCOi_CONSTEXPR DCO_T lowest() { return numeric_limits<typename DCO_T::VALUE_TYPE>::lowest(); }
#endif
    };
}
#endif

