#ifndef _JCALC2_H
#define	_JCALC2_H

#include <fock.h>


namespace votca { namespace moo {

class JCalc
{
public:

    JCalc() { };
   ~JCalc() { };

    void                Initialize() { cout << "MOO"; }
    void                CalcJ() { cout << "MOO"; }
    void                WriteDiPro() { cout << "MOO"; }

private:


    




    struct JCalcData
    {
        mol_and_orb _mol1;
        mol_and_orb _mol2;

        orb         _orb1;
        orb         _orb2;

        fock        _fock;
    };



};




}}

#endif
