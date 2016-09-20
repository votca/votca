/* 
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _JCALC2_H
#define	_JCALC2_H

#include <fock.h>


namespace votca { namespace xtp {

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
