#ifndef FILE_MOOGLOB
#define FILE_MOOGLOB


static double RA=0.529189379; // the bohr radius
const double PI = 4.0 *atan(1.);
const double hbar = 6.582e-16;
const double hbar_si = 1.05457148e-34;
const double Ang = 1E-10;
const double Hartree =27.21; //convert H to eV
/* this variable should sit in global.h*/
        
template<typename T> void clearListList(T&);

#endif //FILE_GLOB
