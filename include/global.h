#ifndef FILE_GLOB
#define FILE_GLOB

#include <iostream>

using namespace std;

static double RA=0.529189379; // the bohr radius
const double PI = 4.0 *atan(1.);
const double hbar = 6.582e-16;
const double hbar_si = 1.05457148e-34;
const double Ang = 1E-10;
const double Hartree =27.21; //convert H to eV
/* this variable should sit in global.h*/
        

template<typename T>
inline void clearListList(T&obj)
{
    typename T::iterator it = obj.begin();
    for ( ; it != obj.end() ; ++it){
        it->clear();
    }
    obj.clear();
}

template<typename T>
inline void safe_delete(T& obj)
{
        typename T::iterator it =obj.begin();
        for( ; it != obj.end() ; ++it ){
//                cout << "call safe_delete" <<endl;
                delete *it;
 //               cout << "called safe_delete" <<endl;
        }
  //      cout << "About to clear safedelete" <<endl;
        obj.clear();
  //      cout << "Cleared safedelete" <<endl;
}


#endif //FILE_GLOB
