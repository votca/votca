/* 
 * File:   units.h
 * Author: vehoff
 *
 * Created on January 21, 2010, 12:07 PM
 */

#ifndef _UNITS_H
#define	_UNITS_H

/// defines all constants necessary to convert from MD to QM units and back
static double RA=0.529189379; /// the Bohr radius
static double hbar_eV=6.58211899e-16; /// Planck constant in [eV*s]
static double Ang = 1E-10; /// Angstroem in [m]
static double Hartree =27.21; /// Hartree in [eV]

/// units allowed in conversion
enum {
    nm=0,
    bohr,
    eV,
    hartree
};

/// template for compile time assertion required for unit conversion class
template <bool t>
struct ctassert {
  enum { N = 1 - 2 * int(!t) };
  /// 1 if t is true, -1 if t is false.
  static char A[N];
  /// ***NOTICE: If this error occurs you tried to use an undefined conversion ***
};

/// class for unit conversion
template<int,int>
class unit{
public:
    template <typename T>
    static T to(T d){
        ctassert<false> foo;
        return T();
    }
};

/// macro to allow easy conversion definitions
#define DEFINE_CONVERSION(a, b, c) \
template <> \
class unit <a,b> {\
public: \
    template <typename T>\
    static T to(T value){\
        return c;\
    }\
};

DEFINE_CONVERSION(nm, bohr, value*10./RA);

#endif /* _UNITS_H */