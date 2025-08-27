/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#ifndef DCO_COMPLEX_HPP
#define DCO_COMPLEX_HPP

//** internal complex type (used for specialization of std::complex later)
namespace dco { namespace helper {
template<typename T>
  class complex_t {
private:
  T _real;
  T _imag;
public:

  DCO_STRONG_INLINE complex_t(T const& real = 0.0, T const& imag = 0.0)
    :_real(real), _imag(imag) { }

  template <typename U>
    complex_t(std::complex<U> const& z) {
    _real = z.real();
    _imag = z.imag();
  }

  DCO_STRONG_INLINE complex_t operator-() const {
    return complex_t(-_real, -_imag);
  }

  DCO_STRONG_INLINE complex_t& operator+=(complex_t<T> const& other) {
    _real += other._real;
    _imag += other._imag;
    return *this;
  }

  DCO_STRONG_INLINE complex_t& operator-=(complex_t const& other) {
    return *this += -other;
  }

  DCO_STRONG_INLINE complex_t& operator*=(complex_t const& other) {
    T const real = _real * other._real - _imag * other._imag;
    T const imag = _real * other._imag + _imag * other._real;
    _real = real;
    _imag = imag;
    return *this;
  }

  DCO_STRONG_INLINE complex_t& operator/=(complex_t const& other) {
    complex_t x1 = (*this) * complex_t(other._real, -other._imag);
    T x2 = other._real * other._real + other._imag * other._imag;
    _real = x1._real / x2;
    _imag = x1._imag / x2;
    return *this;
  }

  DCO_STRONG_INLINE complex_t operator+(complex_t const& other) const {
    complex_t ret = *this;
    return ret += other;
  }
  DCO_STRONG_INLINE complex_t operator-(complex_t const& other) const {
    complex_t ret = *this;
    return ret -= other;
  }
  DCO_STRONG_INLINE complex_t operator*(complex_t const& other) const {
    complex_t ret = *this;
    return ret *= other;
  }
  DCO_STRONG_INLINE complex_t operator/(complex_t const& other) const {
    complex_t ret = *this;
    return ret /= other;
  }

  template <typename U>
    complex_t& operator=(U const& rhs) {
    _real = rhs;
    _imag = T();
    return *this;
  }

  template <typename U>
    complex_t& operator=(std::complex<U> const& rhs) {
    _real = rhs.real();
    _imag = rhs.imag();
    return *this;
  }
  
  DCO_STRONG_INLINE T real() const { return _real; }
  DCO_STRONG_INLINE T imag() const { return _imag; }

  DCO_STRONG_INLINE T& real() { return _real; }
  DCO_STRONG_INLINE T& imag() { return _imag; }

  DCO_STRONG_INLINE void real(T const& real) { _real = real; }
  DCO_STRONG_INLINE void imag(T const& imag) { _imag = imag; }
};

}}

//** specialization of std::complex for dco_type
#ifndef DCO_USE_STD_COMPLEX
namespace std {

template<class DCO_TAPE_REAL, typename DATA_T>
  struct complex<dco::internal::active_type<DCO_TAPE_REAL, DATA_T> >
  : public dco::helper::complex_t<dco::internal::active_type<DCO_TAPE_REAL, DATA_T> >
{
  typedef dco::internal::active_type<DCO_TAPE_REAL, DATA_T> active_t;
  typedef dco::helper::complex_t<active_t >                 dco_complex_t;

  complex(active_t const& real = active_t(),
          active_t const& imag = active_t())
    : dco_complex_t(real, imag) { }
    
  template<typename U>
    complex(complex<U> const& z)
      : dco_complex_t(z) { }

  complex& operator=(active_t const& rhs) {
    static_cast<dco_complex_t&>(*this) = rhs;
    return *this;
  }

  template <typename U>
    complex& operator=(complex<U> const& rhs) {
    static_cast<dco_complex_t&>(*this) = rhs;
    return *this;
  }
};

#define DCO_TEMPLATE_DEF template<class DCO_TAPE_REAL, typename DATA_T> DCO_STRONG_INLINE
#define DCO_TYPE dco::internal::active_type<DCO_TAPE_REAL, DATA_T> 

//** forward declares
DCO_TEMPLATE_DEF DCO_TYPE          abs (complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF DCO_TYPE          arg (complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> log (complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> sin (complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> cos (complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> exp (complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> cosh(complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> tanh(complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> sqrt(complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> pow (complex<DCO_TYPE> const&,
                                                        complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> pow (DCO_TYPE          const&,
                                                        complex<DCO_TYPE> const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> pow (complex<DCO_TYPE> const&,
                                                        DCO_TYPE          const&);
DCO_TEMPLATE_DEF complex<DCO_TYPE> polar(DCO_TYPE const&, const DCO_TYPE& = 0);

//** implementation
DCO_TEMPLATE_DEF
  complex<DCO_TYPE>
  operator+(complex<DCO_TYPE> const& x1, complex<DCO_TYPE> const& x2) {
  return complex<DCO_TYPE>(x1.real() + x2.real(), x1.imag() + x2.imag());
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE>
  operator+(DCO_TYPE const& x1, complex<DCO_TYPE> const& x2) {
  return complex<DCO_TYPE>(x1 + x2.real(), x2.imag());
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE>
  operator+(complex<DCO_TYPE> const& x1, DCO_TYPE const& x2) {
  return complex<DCO_TYPE>(x1.real() + x2, x1.imag());
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE>
  operator-(complex<DCO_TYPE> const& x1, complex<DCO_TYPE> const& x2) {
  return complex<DCO_TYPE>(x1.real() - x2.real(), x1.imag() - x2.imag());
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE>
  operator-(complex<DCO_TYPE> const& x1, DCO_TYPE const& x2) {
  return complex<DCO_TYPE>(x1.real() - x2, x1.imag());
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE>
  operator-(DCO_TYPE const& x1, complex<DCO_TYPE> const& x2) {
  return complex<DCO_TYPE>(x1 - x2.real(), -x2.imag());
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> operator*(complex<DCO_TYPE> const& x1, complex<DCO_TYPE> const& x2) {
  DCO_TYPE real = x1.real() * x2.real() - x1.imag() * x2.imag();
  DCO_TYPE imag = x1.real() * x2.imag() + x1.imag() * x2.real();
  return complex<DCO_TYPE>(real, imag);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> operator*(complex<DCO_TYPE> const& x1, DCO_TYPE const& x2) {
  DCO_TYPE real = x1.real() * x2;
  DCO_TYPE imag = x1.imag() * x2;
  return complex<DCO_TYPE>(real, imag);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> operator*(DCO_TYPE const& x1, complex<DCO_TYPE> const& x2) {
  DCO_TYPE real = x1 * x2.real();
  DCO_TYPE imag = x1 * x2.imag();
  return complex<DCO_TYPE>(real, imag);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> operator/(complex<DCO_TYPE> const& x1, complex<DCO_TYPE> const& x2) {
  complex<DCO_TYPE> a = x1 * complex<DCO_TYPE>(x2.real(), -x2.imag());
  DCO_TYPE b = x2.real() * x2.real() + x2.imag() * x2.imag();
  return complex<DCO_TYPE>(a.real() / b, a.imag() / b);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> operator/(complex<DCO_TYPE> const& x1, DCO_TYPE const& x2) {
  complex<DCO_TYPE> a = x1 * x2;
  DCO_TYPE b = x2 * x2;
  return complex<DCO_TYPE>(a.real() / b, a.imag() / b);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> operator/(DCO_TYPE const& x1, complex<DCO_TYPE> const& x2) {
  complex<DCO_TYPE> a = x1 * complex<DCO_TYPE>(x2.real(), -x2.imag());
  DCO_TYPE b = x2.real() * x2.real() + x2.imag() * x2.imag();
  return complex<DCO_TYPE>(a.real() / b, a.imag() / b);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> operator/(complex<DCO_TYPE> const& x1, const DCO_TAPE_REAL& x2) {
  return x1 / complex<DCO_TYPE>(x2, 0);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> polar(DCO_TYPE const& rho, DCO_TYPE const& theta) {
  return complex<DCO_TYPE>(rho * cos(theta), rho * sin(theta));
}

DCO_TEMPLATE_DEF
  DCO_TYPE abs(complex<DCO_TYPE> const& x) {
  return sqrt( pow(x.real(),2) + pow(x.imag(),2) );
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> cosh(complex<DCO_TYPE> const& x) {
  DCO_TYPE const& x1 = x.real();
  DCO_TYPE const& x2 = x.imag();
  return complex<DCO_TYPE>(cosh(x1) * cos(x2), sinh(x1) * sin(x2));
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> cos(complex<DCO_TYPE> const& x) {
  DCO_TYPE const& x1 = x.real();
  DCO_TYPE const& x2 = x.imag();
  return complex<DCO_TYPE>(cos(x1) * cosh(x2), -sin(x1) * sinh(x2));
}


DCO_TEMPLATE_DEF
  complex<DCO_TYPE> sqrt(complex<DCO_TYPE> const& x) {
  DCO_TYPE const& x1 = x.real();
  DCO_TYPE const& x2 = x.imag();

  if (x1 == DCO_TAPE_REAL()) {
      DCO_TYPE tmp = sqrt(abs(x2) / 2);
      return complex<DCO_TYPE>(tmp, x2 < 0.0 ? DCO_TYPE(-tmp) : tmp);
  } else {
    DCO_TYPE tmp1 = sqrt(2 * (abs(x) + abs(x1)));
    DCO_TYPE tmp2 = tmp1 / 2;
    return x1 > DCO_TAPE_REAL()
      ? complex<DCO_TYPE>(tmp2, DCO_TYPE(x2 / tmp1))
      : complex<DCO_TYPE>(abs(x2) / tmp1, x2 < DCO_TAPE_REAL() ? DCO_TYPE(-tmp2) : tmp2);
  }
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> pow(complex<DCO_TYPE> const& x1,
                                complex<DCO_TYPE> const& x2) {
  return x1 == DCO_TYPE() ? complex<DCO_TYPE>() : exp(x2 * log(x1));
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> pow(DCO_TYPE          const& x1,
                                complex<DCO_TYPE> const& x2) {

  return x1 > DCO_TYPE() ? polar(static_cast<const DCO_TYPE>(pow(x1, x2.real())),
                                 static_cast<const DCO_TYPE>(x2.imag() * log(x1)))
    : pow(complex<DCO_TYPE>(x1), x2);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> pow(complex<DCO_TYPE> const& x1,
                                DCO_TYPE  const& x2) {

  if (x1.imag() == DCO_TAPE_REAL() && x1.real() > DCO_TAPE_REAL())
    return static_cast<DCO_TYPE>(pow(x1.real(), x2));
  
  complex<DCO_TYPE> t = log(x1);
  return polar(static_cast<const DCO_TYPE>(exp(x2 * t.real())),
               static_cast<const DCO_TYPE>(x2 * t.imag()));
}


DCO_TEMPLATE_DEF
  complex<DCO_TYPE> exp(complex<DCO_TYPE> const& x) {
  return polar(static_cast<const DCO_TYPE>(exp(x.real())), x.imag());
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> log(complex<DCO_TYPE> const& x) {
  return complex<DCO_TYPE>(log(abs(x)), arg(x));
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> log10(complex<DCO_TYPE> const& x) {
  return complex<DCO_TYPE>(log(x) / log(DCO_TAPE_REAL(10.0)));
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> sin(complex<DCO_TYPE> const& x) {
  DCO_TYPE const& x1 = x.real();
  DCO_TYPE const& x2 = x.imag();
  return complex<DCO_TYPE>(sin(x1) * cosh(x2), cos(x1) * sinh(x2));
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> sinh(complex<DCO_TYPE> const& x) {
  DCO_TYPE const&  x1 = x.real();
  DCO_TYPE const&  x2 = x.imag();
  return complex<DCO_TYPE>(sinh(x1) * cos(x2), cosh(x1) * sin(x2));
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> tan(complex<DCO_TYPE> const& x) {
  return sin(x) / cos(x);
}

DCO_TEMPLATE_DEF
  complex<DCO_TYPE> tanh(complex<DCO_TYPE> const& x) {
  return sinh(x) / cosh(x);
}

DCO_TEMPLATE_DEF
  DCO_TYPE arg(complex<DCO_TYPE> const& x) {
  return  atan2(x.imag(), x.real());
}

#undef DCO_TEMPLATE_DEF
#undef DCO_TYPE

} // end namespace std
#endif //DCO_USE_STD_COMPLEX
#endif // DCO_COMPLEX_HPP


