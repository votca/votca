/* 
 * Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _atom_table_H
#define _atom_table_H

#include <vector>
#include <string>
#include <unordered_map>

namespace votca { namespace csg {

/**
 * Brief information about atom table
 *
 * The atom table class is responsible for stroring relevant information 
 * from the periodic table.
 *
 */

class AtomTable{

    public:
        /**
         * Constructor
         */
        AtomTable(void);
        /* 
         * Destructor 
         */
        ~AtomTable() {}

        /* *
         * Check to see if the symbol is a valid symbol. The symbols follow the
         * standard seen in the periodic table, H - Hydrogen, He - Helium etc.
         *
         * Returns:
         * 
         * true  - if symbol is real
         * 
         * false - if symbol is not recognized
         */
        bool checkSymbol(std::string sym);
        /**
         * Will return the symbols associated with the halogen atoms in a vector<string>.
         */
        std::vector<std::string> getHalogens(void);
        /**
         *  Will return the symbols associated with the noble gasses in a vector<string>.
         */
        std::vector<std::string> getNoble(void);

        /**
         * Pass the symbol in and returns the mass in a.u.
         *
         * If the symbol is not recognized will exit with an error message
         *
         * If the mass is not known it will also exit with an error message
         */
        double getMass(std::string sym);
        /**
         * Pass the symbol in and returns the atomic number.
         *
         * If the symbol is not recognized it will print an error message and exit.
         */
        int getAtomicNumber(std::string sym);

    protected:

        struct Atom {
            double mass;
            int atomicNumber;
        };

        std::unordered_map<std::string, Atom> atomMap;
        /**
         * Return the index associated with the tables for the chosen element.
         */
        int getIndex(std::string sym);

        /**
         * Stores the symbols for all the atoms in the periodic table.
         */
    std::vector<std::string> Atoms = {"Ac" ,"Ag" ,"Al" ,"Am" ,\
  "Ar" ,"As" ,"At" ,"Au" ,"B"  ,"Ba" ,"Be" ,"Bh" ,"Bi" ,"Bk" ,\
  "Br" ,"C"  ,"Ca" ,"Cd" ,"Ce" ,"Cf" ,"Cl" ,"Cm" ,"Cn" ,"Co" ,\
  "Cr" ,"Cs" ,"Cu" ,"Db" ,"Ds" ,"Dy" ,"Er" ,"Es" ,"Eu" ,"F"  ,\
  "Fe" ,"Fl" ,"Fm" ,"Fr" ,"Ga" ,"Gd" ,"Ge" ,"H"  ,"He" ,"Hf" ,\
  "Hg" ,"Ho" ,"Hs" ,"I"  ,"In" ,"Ir" ,"K"  ,"Kr" ,"La" ,"Li" ,\
  "Lr" ,"Lu" ,"Lv" ,"Md" ,"Mg" ,"Mn" ,"Mo" ,"Mt" ,"N"  ,"Na" ,\
  "Nb" ,"Nd" ,"Ne" ,"Ni" ,"No" ,"Np" ,"O"  ,"Os" ,"P"  ,"Pa" ,\
  "Pb" ,"Pd" ,"Pm" ,"Po" ,"Pr" ,"Pt" ,"Pu" ,"Ra" ,"Rb" ,"Re" ,\
  "Rf" ,"Rg" ,"Rh" ,"Rn" ,"Ru" ,"S"  ,"Sb" ,"Sc" ,"Se" ,"Sg" ,\
  "Si" ,"Sm" ,"Sn" ,"Sr" ,"Ta" ,"Tb" ,"Tc" ,"Te" ,"Th" ,"Ti" ,\
  "Tl" ,"Tm" ,"U"  ,"Uub","Uuh","Uuo","Uup","Uuq","Uus","Uut",\
  "V"  ,"W"  ,"Xe" ,"Y"  ,"Yb" ,"Zn" ,"Zr" };

    /**
     * Stores the atomic numbers for all the elements in the periodic table.
     */

    std::vector<int> AtomicNumber  = { 89  , 47  , 13  , 95  ,\
   18  , 33  , 85  , 79  ,  5  , 56  ,  4  , 107 , 83  , 97  ,\
   35  ,  6  , 20  , 48  , 58  , 98  , 17  , 96  , 112 , 27  ,\
   24  , 55  , 29  , 105 , 110 , 66  , 68  , 99  , 63  ,  9  ,\
   26  , 114 , 100 , 87  , 31  , 64  , 32  ,  1  ,  2  , 72  ,\
   80  , 67  , 108 , 53  , 49  , 77  , 19  , 36  , 57  ,  3  ,\
   103 , 71  , 116 , 101 , 12  , 25  , 42  , 109 ,  7  , 11  ,\
   41  , 60  , 10  , 28  , 102 , 93  ,  8  , 76  , 15  , 91  ,\
   82  , 46  , 61  , 84  , 59  , 78  , 94  , 88  , 37  , 75  ,\
   104 , 111 , 45  , 86  , 44  , 16  , 51  , 21  , 34  , 106 ,\
   14  , 62  , 50  , 38  , 73  , 65  , 43  , 52  , 90  , 22  ,\
   81  , 69  , 92  , 112 , 116 , 118 , 115 , 114 , 117 , 113 ,\
   23  , 74  , 54  , 39  , 70  , 30  , 40  };

    /**
     * Stores the masses for all the elements in the periodic table.
     *
     * A value of -1 indicates the mass is unknown for the element.
     */
    std::vector<double> MassNumber =         {227.028  ,107.8682  , 26.981539,241.056829,\
    39.948  , 74.92159 ,209.987148,169.96654 , 10.811  ,137.327   ,  9.012182,272.1380  ,208.98037 ,247.070307  ,
    79.904  , 12.011   , 40.078   ,112.411   ,140.115  ,249.074853, 35.4527  ,243.061389,285.1741  , 58.93320   ,
    51.9961 ,132.90543 , 63.546   ,268.1255  ,281.1621 ,162.50    ,167.26    ,252.08298 ,151.965   , 18.9984032 ,
    55.847  ,289.1873  ,257.095105,223.019736, 69.723  ,157.25    , 72.61    ,  1.00794 ,  4.002602,178.49      ,
   200.59   ,164.93032 ,270.13    ,126.90447 ,114.82   ,192.22    , 39.0983  , 83.80    ,138.9055  ,  6.941     ,
   262.1096 ,174.967   , -1       ,258       , 24.3050 , 54.93805 , 95.94    ,276.1512  , 14.00674 , 22.98976928,
    92.90638,144.24    , 20.1797  , 58.6934  ,259.1010 ,237.048   , 15.9994  ,190.2     , 30.973762,231.0359    ,
   207.2    ,106.42    ,144.912749,209       ,140.90765,195.08    ,238.04956 ,266.025   , 85.4678  ,186.207     ,
   265.1167 ,280.1645  ,102.90550 ,222       ,101.07   , 32.066   ,121.760   , 44.955910, 78.96    ,271.1335    ,
    28.0855 ,150.36    ,118.710   , 87.62    ,180.9479 ,158.92534 , 96.906365,127.60    ,232.0381  , 47.88      ,
   204.3833 ,168.93421 ,238.0289  ,285.1741  , -1      , -1       ,288.1925  ,289.1873  ,292.208   ,284.1781    ,
    50.9415 ,183.85    ,131.29    , 88.90585 ,173.04   , 65.39    , 91.224   };  
     
};
}}
#endif
