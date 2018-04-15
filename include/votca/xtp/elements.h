/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef __VOTCA_XTP_ELEMENTS_H
#define	__VOTCA_XTP_ELEMENTS_H

#include <string>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>



namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;



/**
    \brief information about an element
 
    The Atom class stores atom id, name, type, mass, charge, residue number
    
*/
class Elements
{
public:   

  

    Elements() { FillMaps(); };
   ~Elements() { };

    const double     &getVdWChelpG(std::string name ) const { return _VdWChelpG.at(name); }
    const double     &getVdWMK(std::string name ) const { return _VdWMK.at(name); }
    const double     &getNucCrgECP(std::string name) const {return _NucCrgECP.at(name); }
    const double     &getNucCrg(std::string name) const {return _NucCrg.at(name); }
    const int        &getEleNum(std::string name) const {return _EleNum.at(name); }
    const double     &getMass(std::string name) const {return _Mass.at(name); }

    const std::string     &getEleName(int elenum) const {return _EleName.at(elenum); }
    const std::string     &getEleShort(std::string elefull) const {return _EleShort.at(elefull); }
    const std::string     &getEleFull(std::string eleshort) const {return _EleFull.at(eleshort); }
    
   
private:

    std::map<std::string, double> _VdWChelpG;
    std::map<std::string, double> _VdWMK;
    std::map<std::string, double> _NucCrgECP;
    std::map<std::string, double> _NucCrg;
    std::map<std::string, double> _Mass;
    std::map<std::string, int> _EleNum;
    std::map<int, std::string> _EleName;
    
    std::map< std::string,ub::matrix<int> > _Atomconfig;
    
    std::map<std::string, std::string> _EleShort;
    std::map<std::string, std::string> _EleFull;
    
    inline void FillMaps(){
        
        FillVdWChelpG();
        FillVdWMK();
        FillNucCrgECP();
        FillNucCrg();
        FillEleNum();
        FillEleName();
        FillEleShort();
        FillEleFull();       
        FillMass();
    };
    
    
    
    inline void FillMass(){
    
        // masses of atoms
        _Mass["H"]  = 1.00794; 
        _Mass["He"] = 4.002602;
        _Mass["Li"] = 6.941;
        _Mass["Be"] = 9.012182;
        _Mass["B"]  = 10.811;
        _Mass["C"]  = 12.0107;
        _Mass["N"]  = 14.00674;
        _Mass["O"]  = 15.9994;
        _Mass["F"]  = 18.9984032;
        _Mass["Ne"] = 20.1797;
        _Mass["Na"] = 22.989770;
        _Mass["Mg"] = 24.3050;
        _Mass["Al"] = 26.981538;
        _Mass["Si"] = 28.0855;
        _Mass["P"]  = 30.973761;
        _Mass["S"]  = 32.066;
        _Mass["Cl"] = 35.4527;
        _Mass["Ar"] = 39.948;
        _Mass["Hg"] = 200.59;
        _Mass["Ag"] = 107.8682;
    };

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
        _VdWChelpG["Hg"] = 2.0;
        _VdWChelpG["Ag"] = 1.7;
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
        _NucCrgECP["Ag"] = 19.00;
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
        _NucCrg["Ag"] = 47.00;
    };
    
    
     inline void FillEleNum(){
    
        // Nuclear Charges
        _EleNum["H"]  = 1; 
        _EleNum["He"] = 2;
        _EleNum["Li"] = 3;
        _EleNum["Be"] = 4;
        _EleNum["B"]  = 5;
        _EleNum["C"]  = 6;
        _EleNum["N"]  = 7;
        _EleNum["O"]  = 8;
        _EleNum["F"]  = 9;
        _EleNum["Ne"] = 10;
        _EleNum["Na"] = 11;
        _EleNum["Mg"] = 12;
        _EleNum["Al"] = 13;
        _EleNum["Si"] = 14;
        _EleNum["P"]  = 15;
        _EleNum["S"]  = 16;
        _EleNum["Cl"] = 17;
        _EleNum["Ar"] = 18;
        _EleNum["Ag"] = 47;
    };
    
    
     inline void FillEleName(){
    
        // Nuclear Charges
        _EleName[1]  = "H"; 
        _EleName[2] = "He";
        _EleName[3] = "Li";
        _EleName[4] = "Be";
        _EleName[5]  ="B";
        _EleName[6]  = "C";
        _EleName[7]  = "N";
        _EleName[8]  = "O";
        _EleName[9]  = "F";
        _EleName[10] = "Ne";
        _EleName[11] = "Na";
        _EleName[12] = "Mg";
        _EleName[13] = "Al";
        _EleName[14] = "Si";
        _EleName[15]  = "P";
        _EleName[16]  = "S";
        _EleName[17] = "Cl";
        _EleName[18] = "Ar";
        _EleName[47] = "Ag";
    };
    
    inline void FillEleShort(){
    
        // VdW radii in Angstrom as used in MK Gaussian
        _EleShort["HYDROGEN"]  = "H"; 
        _EleShort["HELIUM"] = "He";
        _EleShort["LITHIUM"] = "Li";
        _EleShort["BERYLLIUM"] = "Be";
        _EleShort["BORON"]  = "B";
        _EleShort["CARBON"]  = "C";
        _EleShort["NITROGEN"]  = "Ni";
        _EleShort["OXYGEN"]  = "O";
        _EleShort["FLUORINE"]  = "Fl";
        _EleShort["NEON"] = "Ne";
        _EleShort["SODIUM"] = "Na";
        _EleShort["MAGNESIUM"] = "Mg";
        _EleShort["ALUMINUM"] = "Al";
        _EleShort["SILICON"] = "Si";
        _EleShort["PHOSPHORUS"]  = "Ph";
        _EleShort["SULFUR"]  = "S";
        _EleShort["CLORINE"] = "Cl";
        _EleShort["ARGON"] = "Ar";
        _EleShort["SILVER"] = "Ag";        
    };
    
 inline void FillEleFull(){
    
        // VdW radii in Angstrom as used in MK Gaussian
        _EleFull["H"]  = "HYDROGEN"; 
        _EleFull["He"] = "HELIUM";
        _EleFull["Li"] = "LITHIUM";
        _EleFull["Be"] = "BERYLLIUM";
        _EleFull["B"]  = "BORON";
        _EleFull["C"]  = "CARBON";
        _EleFull["N"]  = "NITROGEN";
        _EleFull["O"]  = "OXYGEN";
        _EleFull["F"]  = "FLUORINE";
        _EleFull["Ne"] = "NEON";
        _EleFull["Na"] = "SODIUM";
        _EleFull["Mg"] = "MAGNESIUM";
        _EleFull["Al"] = "ALUMINUM";
        _EleFull["Si"] = "SILICON";
        _EleFull["P"]  = "PHOSPHORUS";
        _EleFull["S"]  = "SULFUR";
        _EleFull["Cl"] = "CHLORINE";
        _EleFull["Ar"] = "ARGON";
        _EleFull["Ag"] = "SILVER";        
    };
   

    
    
    inline void FillVdWMK(){
    
        // VdW radii in Angstrom as used in MK Gaussian
        _VdWMK["H"]  = 1.2; 
        _VdWMK["He"] = 1.2;
        _VdWMK["Li"] = 1.37;
        _VdWMK["Be"] = 1.45;
        _VdWMK["B"]  = 1.45;
        _VdWMK["C"]  = 1.5;
        _VdWMK["N"]  = 1.5;
        _VdWMK["O"]  = 1.4;
        _VdWMK["F"]  = 1.35;
        _VdWMK["Ne"] = 1.3;
        _VdWMK["Na"] = 1.57;
        _VdWMK["Mg"] = 1.36;
        _VdWMK["Al"] = 1.24;
        _VdWMK["Si"] = 1.17;
        _VdWMK["P"]  = 1.8;
        _VdWMK["S"]  = 1.75;
        _VdWMK["Cl"] = 1.7;
        // _VdWMK["Ar"] = 2.0;
        _VdWMK["Ag"] = 2.0;
    };
    
    
    
   
};
}}

#endif	/* __VOTCA_XTP_ATOM_H */

