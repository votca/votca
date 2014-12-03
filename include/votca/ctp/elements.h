/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef __VOTCA_CTP_ELEMENTS_H
#define	__VOTCA_CTP_ELEMENTS_H

#include <string>
#include <assert.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>
#include <votca/tools/property.h>
#include <votca/tools/matrix.h>
#include <map>



namespace votca { namespace ctp {
using namespace votca::tools;

using namespace std;


/**
    \brief information about an element
 
    The Atom class stores atom id, name, type, mass, charge, residue number
    
*/
class Elements
{
public:   

  

    Elements() { FillMaps(); };
   ~Elements() { };

    const double     &getVdWChelpG(string name ) const { return _VdWChelpG.at(name); }
    const double     &getVdWMK(string name ) const { return _VdWMK.at(name); }
    const double     &getNucCrgECP(string name) const {return _NucCrgECP.at(name); }
    const double     &getNucCrg(string name) const {return _NucCrg.at(name); }

private:

    std::map<std::string, double> _VdWChelpG;
    std::map<std::string, double> _VdWMK;
    std::map<std::string, double> _NucCrgECP;
    std::map<std::string, double> _NucCrg;
    
    inline void FillMaps(){
        
        FillVdWChelpG();
        FillVdWMK();
        FillNucCrgECP();
        FillNucCrg();
        
    }

    inline void FillVdWChelpG(){
    
        // VdW radii in Angstrom as used in CHELPG paper [Journal of Computational Chemistry 11, 361, 1990]and Gaussian
        _VdWChelpG["H"]  = 1.45; 
        _VdWChelpG["He"] = 1.45;
        _VdWChelpG["Li"] = 1.5;
        _VdWChelpG["Be"] = 1.5;
        _VdWChelpG["B"]  = 1.5;
        _VdWChelpG["C"]  = 1.5;
        _VdWChelpG["N"]  = 1.7;
        _VdWChelpG["O"]  = 1.7;
        _VdWChelpG["F"]  = 1.7;
        _VdWChelpG["Ne"] = 1.7;
        _VdWChelpG["Na"] = 2.0;
        _VdWChelpG["Mg"] = 2.0;
        _VdWChelpG["Al"] = 2.0;
        _VdWChelpG["Si"] = 2.0;
        _VdWChelpG["P"]  = 2.0;
        _VdWChelpG["S"]  = 2.0;
        _VdWChelpG["Cl"] = 2.0;
        _VdWChelpG["Ar"] = 2.0;
    };


    inline void FillNucCrgECP(){
    
        // Effective Core Potential Nuclear Charges
        _NucCrgECP["H"]  = 1.00; 
        _NucCrgECP["He"] = 2.00;
        _NucCrgECP["Li"] = 1.00;
        _NucCrgECP["Be"] = 2.00;
        _NucCrgECP["B"]  = 3.00;
        _NucCrgECP["C"]  = 4.00;
        _NucCrgECP["N"]  = 5.00;
        _NucCrgECP["O"]  = 6.00;
        _NucCrgECP["F"]  = 7.00;
        _NucCrgECP["Ne"] = 8.00;
        _NucCrgECP["Na"] = 1.00;
        _NucCrgECP["Mg"] = 2.00;
        _NucCrgECP["Al"] = 3.00;
        _NucCrgECP["Si"] = 4.00;
        _NucCrgECP["P"]  = 5.00;
        _NucCrgECP["S"]  = 6.00;
        _NucCrgECP["Cl"] = 7.00;
        _NucCrgECP["Ar"] = 8.00;
    };
    
      inline void FillNucCrg(){
    
        // Nuclear Charges
        _NucCrg["H"]  = 1.00; 
        _NucCrg["He"] = 2.00;
        _NucCrg["Li"] = 3.00;
        _NucCrg["Be"] = 4.00;
        _NucCrg["B"]  = 5.00;
        _NucCrg["C"]  = 6.00;
        _NucCrg["N"]  = 7.00;
        _NucCrg["O"]  = 8.00;
        _NucCrg["F"]  = 9.00;
        _NucCrg["Ne"] = 10.00;
        _NucCrg["Na"] = 11.00;
        _NucCrg["Mg"] = 12.00;
        _NucCrg["Al"] = 13.00;
        _NucCrg["Si"] = 14.00;
        _NucCrg["P"]  = 15.00;
        _NucCrg["S"]  = 16.00;
        _NucCrg["Cl"] = 17.00;
        _NucCrg["Ar"] = 18.00;
    };
    
    
    inline void FillVdWMK(){
    
        // VdW radii in Angstrom as used in MK Gaussian
        _VdWChelpG["H"]  = 1.2; 
        _VdWChelpG["He"] = 1.2;
        _VdWChelpG["Li"] = 1.37;
        _VdWChelpG["Be"] = 1.45;
        _VdWChelpG["B"]  = 1.45;
        _VdWChelpG["C"]  = 1.5;
        _VdWChelpG["N"]  = 1.5;
        _VdWChelpG["O"]  = 1.4;
        _VdWChelpG["F"]  = 1.35;
        _VdWChelpG["Ne"] = 1.3;
        _VdWChelpG["Na"] = 1.57;
        _VdWChelpG["Mg"] = 1.36;
        _VdWChelpG["Al"] = 1.24;
        _VdWChelpG["Si"] = 1.17;
        _VdWChelpG["P"]  = 1.8;
        _VdWChelpG["S"]  = 1.75;
        _VdWChelpG["Cl"] = 1.7;
        // _VdWChelpG["Ar"] = 2.0;
    };
};




}}

#endif	/* __VOTCA_CTP_ATOM_H */

