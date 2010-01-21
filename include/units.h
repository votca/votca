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

enum {
    u_nm=0,
    u_bohr,
    u_eV,
    u_hartree
};

template <bool t>
struct ctassert {
  enum { N = 1 - 2 * int(!t) };
  /// 1 if t is true, -1 if t is false.
  static char A[N];
  /// ***NOTICE: If this error occurs you tried to use an undefined conversion ***
};

template<int, int>
double convert(const double &value)
{
   // compile time assetion
    ctassert<false> foo;
}
/// convert from nm to Bohr
template<>
inline double convert<u_nm, u_bohr>(const double &value)
{
    return value*RA/10.;
}
/// convert from Bohr to nm
template<>
inline double convert<u_bohr, u_nm>(const double &value)
{
    return value*10./RA;
}
/// convert from eV to Hartree
template<>
inline double convert<u_eV, u_hartree>(const double &value)
{
    return value/Hartree;
}
/// convert from Hartree to eV
template<>
inline double convert<u_hartree, u_eV>(const double &value)
{
    return value*Hartree;
}
#endif	/* _UNITS_H */

