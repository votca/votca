/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef VOTCA_XTP_EULER_MACLAURIN_H
#define	VOTCA_XTP_EULER_MACLAURIN_H


#include <votca/xtp/basisset.h>
#include <votca/xtp/qmatom.h>
#include <votca/tools/constants.h>
#include <votca/xtp/grid_containers.h>

namespace votca { namespace xtp {

    
        class EulerMaclaurinGrid {
        public: 
            
            EulerMaclaurinGrid() { FillGrids(); };
            
            std::map<std::string, GridContainers::radial_grid> CalculateAtomicRadialGrids(const AOBasis& aobasis,
                                                        const QMMolecule& atoms , const std::string& type);
            
            std::vector<double> CalculatePruningIntervals( const std::string& element );
            
        private:
            
            struct min_exp {
                double alpha;
                int l;
                double range;
            };
       
            std::map<std::string,min_exp> _element_ranges;
            std::map<std::string,int> _pruning_set;
            
            std::map<std::string,double> _BraggSlaterRadii;
            
            int getGridParameters(const std::string& element, const std::string& type);
            
            double DetermineCutoff( double alpha, int l, double eps );
            double CalcResidual( double alpha, int l, double cutoff);
            double RadialIntegral(double alpha, int l, double cutoff);
            
            void CalculateRadialCutoffs(const AOBasis& aobasis,const QMMolecule& atoms , const std::string& gridtype );
            void RefineElementRangeMap(const AOBasis& aobasis, const QMMolecule& atoms , double eps);

            void FillElementRangeMap(const AOBasis& aobasis, const QMMolecule& atoms , double eps);

            GridContainers::radial_grid CalculateRadialGridforAtom(const std::string& type,const std::pair<std::string,min_exp>& element);
            
            std::map<std::string, int>    MediumGrid;
            std::map<std::string, int>    CoarseGrid;
            std::map<std::string, int>    XcoarseGrid;
            std::map<std::string, int>    FineGrid;
            std::map<std::string, int>    XfineGrid;
            std::map<std::string, double> Accuracy;
            
            inline void FillGrids(){  
                FillBraggSlaterRadii();
                FillPruningSet();
                FillAccuracy();
                FillMediumGrid();
                FillCoarseGrid();
                FillXcoarseGrid();
                FillFineGrid();
                FillXfineGrid();
            }

            
            inline void FillBraggSlaterRadii(){
                const double ang2bohr=votca::tools::conv::ang2bohr; 
                // 
                _BraggSlaterRadii["H"] = 0.35 * ang2bohr ;
                _BraggSlaterRadii["He"] = 0.35 * ang2bohr ;

                // row of period system  for 1st row elements taken from NWChem
                _BraggSlaterRadii["Li"] = 1.45 * ang2bohr ;
                _BraggSlaterRadii["Be"] = 1.05* ang2bohr ;
                _BraggSlaterRadii["B"] = 0.85* ang2bohr ;
                _BraggSlaterRadii["C"] = 0.70* ang2bohr ;
                _BraggSlaterRadii["N"] = 0.65* ang2bohr ;
                _BraggSlaterRadii["O"] = 0.60* ang2bohr ;
                _BraggSlaterRadii["F"] = 0.50* ang2bohr ;
                _BraggSlaterRadii["Ne"] = 0.50* ang2bohr ;

                // row of period system  for 2nd row elements taken from NWChem
                _BraggSlaterRadii["Na"] = 1.80* ang2bohr ;
                _BraggSlaterRadii["Mg"] = 1.5* ang2bohr ;
                _BraggSlaterRadii["Al"] = 1.25* ang2bohr ;
                _BraggSlaterRadii["Si"] = 1.1* ang2bohr ;
                _BraggSlaterRadii["P"] = 1.0* ang2bohr ;
                _BraggSlaterRadii["S"] = 1.0* ang2bohr ;
                _BraggSlaterRadii["Cl"] = 1.0* ang2bohr ;
                _BraggSlaterRadii["Ar"] = 1.0* ang2bohr ;

                // row of period system  for 3rd row elements taken from NWChem
                _BraggSlaterRadii["K"] = 2.2* ang2bohr ;
                _BraggSlaterRadii["Ca"] = 1.8* ang2bohr ;
                _BraggSlaterRadii["Sc"] = 1.6* ang2bohr ;
                _BraggSlaterRadii["Ti"] = 1.4* ang2bohr ;
                _BraggSlaterRadii["V"] = 1.35* ang2bohr ;
                _BraggSlaterRadii["Cr"] = 1.4* ang2bohr ;
                _BraggSlaterRadii["Mn"] = 1.4* ang2bohr ;
                _BraggSlaterRadii["Fe"] = 1.4* ang2bohr ;
                _BraggSlaterRadii["Co"] = 1.35* ang2bohr ;
                _BraggSlaterRadii["Ni"] = 1.35* ang2bohr ;
                _BraggSlaterRadii["Cu"] = 1.35* ang2bohr ;
                _BraggSlaterRadii["Zn"] = 1.35* ang2bohr ;
                _BraggSlaterRadii["Ga"] = 1.30* ang2bohr ;
                _BraggSlaterRadii["Ge"] = 1.25* ang2bohr ;
                _BraggSlaterRadii["As"] = 1.15* ang2bohr ;
                _BraggSlaterRadii["Se"] = 1.15* ang2bohr ;
                _BraggSlaterRadii["Br"] = 1.15* ang2bohr ;
                _BraggSlaterRadii["Kr"] = 1.15* ang2bohr ;              
                
                // 4th row (selection)
                _BraggSlaterRadii["Ag"] = 1.60* ang2bohr ;
                
                /* Copied from grid_atom_type_info.F of NWChem
                 
                 * VALUES are in ANGSTROM
                       Data BSrad/0.35,0.35,1.45,1.05,0.85,0.70,0.65,0.60,0.50,0.50,
c                  Na   Mg   Al   Si    P    S   Cl   Ar    K   Ca
     &           1.80,1.50,1.25,1.10,1.00,1.00,1.00,1.00,2.20,1.80,
c                  Sc   Ti    V   Cr   Mn   Fe   Co   Ni   Cu   Zn
     &           1.60,1.40,1.35,1.40,1.40,1.40,1.35,1.35,1.35,1.35,
c                  Ga   Ge   As   Se   Br   Kr   Rb   Sr    Y   Zr
     &           1.30,1.25,1.15,1.15,1.15,1.15,2.35,2.00,1.80,1.55,
c                  Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd   In   Sn
     &           1.45,1.45,1.35,1.30,1.35,1.40,1.60,1.55,1.55,1.45,
c                  Sb   Te    I   Xe   Cs   Ba   La   Ce   Pr   Nd
     &           1.45,1.40,1.40,1.40,2.60,2.15,1.95,1.85,1.85,1.85,
c                  Pm   Sm   Eu   Gd   Tb   Dy   Ho   Er   Tm   Yb
     &           1.85,1.85,1.85,1.80,1.75,1.75,1.75,1.75,1.75,1.75,
c                  Lu   Hf   Ta    W   Re   Os   Ir   Pt   Au   Hg
     &           1.75,1.55,1.45,1.35,1.35,1.30,1.35,1.35,1.35,1.50,
c                  Tl   Pb   Bi   Po   At   Rn   Fr   Ra   Ac   Th
     &           1.90,1.80,1.60,1.90,1.90,1.90,2.60,2.15,1.95,1.80,
c                  Pa    U   Np   Pu   Am   Cm   Bk   Cf   Es   Fm
     &           1.80,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,
c                  Md   No   Lr  Unq  Unp
     &           1.75,1.75,1.75,1.55,1.55/
                 
                 
                 
                 */         
            }
            
            
            
            
            inline void FillPruningSet() {

                // row of period system for H, He (not given in NWChem, assuming same as 1st row)
                _pruning_set["H"] = 1;
                _pruning_set["He"] = 1;

                // row of period system  for 1st row elements taken from NWChem
                _pruning_set["Li"] = 2;
                _pruning_set["Be"] = 2;
                _pruning_set["B"] = 2;
                _pruning_set["C"] = 2;
                _pruning_set["N"] = 2;
                _pruning_set["O"] = 2;
                _pruning_set["F"] = 2;
                _pruning_set["Ne"] = 2;

                // row of period system  for 2nd row elements taken from NWChem
                _pruning_set["Na"] = 3;
                _pruning_set["Mg"] = 3;
                _pruning_set["Al"] = 3;
                _pruning_set["Si"] = 3;
                _pruning_set["P"] = 3;
                _pruning_set["S"] = 3;
                _pruning_set["Cl"] = 3;
                _pruning_set["Ar"] = 3;

                // row of period system  for 3rd row elements taken from NWChem
                _pruning_set["K"] = 3;
                _pruning_set["Ca"] = 3;
                _pruning_set["Sc"] = 3;
                _pruning_set["Ti"] = 3;
                _pruning_set["V"] = 3;
                _pruning_set["Cr"] = 3;
                _pruning_set["Mn"] = 3;
                _pruning_set["Fe"] = 3;
                _pruning_set["Co"] = 3;
                _pruning_set["Ni"] = 3;
                _pruning_set["Cu"] = 3;
                _pruning_set["Zn"] = 3;
                _pruning_set["Ga"] = 3;
                _pruning_set["Ge"] = 3;
                _pruning_set["As"] = 3;
                _pruning_set["Se"] = 3;
                _pruning_set["Br"] = 3;
                _pruning_set["Kr"] = 3;
                
                // row of period system for 4th row elements taken from NWChem
                _pruning_set["Ag"] = 3;
            }
            
            inline void FillAccuracy(){
                Accuracy["xcoarse"] = 1e-4;
                Accuracy["coarse"]  = 1e-5;             
                Accuracy["medium"]  = 1e-6;
                Accuracy["fine"]    = 1e-7;
                Accuracy["xfine"]   = 1e-8;             
            }
            
    inline void FillMediumGrid(){
        
        // order for H, He (not given in NWChem, assuming same as 1st row)
            MediumGrid["H"]  = 49;
            MediumGrid["He"] = 49;
        
        // orders for 1st row elements taken from NWChem
            MediumGrid["Li"] = 49;
            MediumGrid["Be"] = 49;
            MediumGrid["B"] = 49;
            MediumGrid["C"] = 49;
            MediumGrid["N"] = 49;
            MediumGrid["O"] = 49;
            MediumGrid["F"] = 49;
            MediumGrid["Ne"] = 49;

        // orders for 2nd row elements taken from NWChem
            MediumGrid["Na"] = 88;
            MediumGrid["Mg"] = 88;
            MediumGrid["Al"] = 88;
            MediumGrid["Si"] = 88;
            MediumGrid["P"] = 88;
            MediumGrid["S"] = 88;
            MediumGrid["Cl"] = 88;
            MediumGrid["Ar"] = 88;
            
        // orders for 3rd row elements taken from NWChem
            MediumGrid["K"] = 112;
            MediumGrid["Ca"] = 112;
            MediumGrid["Sc"] = 112;
            MediumGrid["Ti"] = 112;
            MediumGrid["V"] = 112;
            MediumGrid["Cr"] = 112;
            MediumGrid["Mn"] = 112;
            MediumGrid["Fe"] = 112;
            MediumGrid["Co"] = 112;
            MediumGrid["Ni"] = 112;
            MediumGrid["Cu"] = 112;
            MediumGrid["Zn"] = 112;
            MediumGrid["Ga"] = 112;
            MediumGrid["Ge"] = 112;
            MediumGrid["As"] = 112;
            MediumGrid["Se"] = 112;
            MediumGrid["Br"] = 112;
            MediumGrid["Kr"] = 112;
            
            // orders for 4th row elements (selected)
            MediumGrid["Ag"] = 123;   
    }            
    
    inline void FillFineGrid(){
        
        // order for H, He (not given in NWChem, assuming same as 1st row)
            FineGrid["H"]  = 70;
            FineGrid["He"] = 70;
        
        // orders for 1st row elements taken from NWChem
            FineGrid["Li"] = 70;
            FineGrid["Be"] = 70;
            FineGrid["B"] = 70;
            FineGrid["C"] = 70;
            FineGrid["N"] = 70;
            FineGrid["O"] = 70;
            FineGrid["F"] = 70;
            FineGrid["Ne"] = 70;

        // orders for 2nd row elements taken from NWChem
            FineGrid["Na"] = 123;
            FineGrid["Mg"] = 123;
            FineGrid["Al"] = 123;
            FineGrid["Si"] = 123;
            FineGrid["P"] = 123;
            FineGrid["S"] = 123;
            FineGrid["Cl"] = 123;
            FineGrid["Ar"] = 123;
            
        // orders for 3rd row elements taken from NWChem
            FineGrid["K"] = 130;
            FineGrid["Ca"] = 130;
            FineGrid["Sc"] = 130;
            FineGrid["Ti"] = 130;
            FineGrid["V"] = 130;
            FineGrid["Cr"] = 130;
            FineGrid["Mn"] = 130;
            FineGrid["Fe"] = 130;
            FineGrid["Co"] = 130;
            FineGrid["Ni"] = 130;
            FineGrid["Cu"] = 130;
            FineGrid["Zn"] = 130;
            FineGrid["Ga"] = 130;
            FineGrid["Ge"] = 130;
            FineGrid["As"] = 130;
            FineGrid["Se"] = 130;
            FineGrid["Br"] = 130;
            FineGrid["Kr"] = 130;
            
            // 4th row (selected)
            FineGrid["Ag"] = 141;
    }            
    
    inline void FillXfineGrid(){
        
        // order for H, He (not given in NWChem, assuming same as 1st row)
            XfineGrid["H"]  = 100;
            XfineGrid["He"] = 100;
        
        // orders for 1st row elements taken from NWChem
            XfineGrid["Li"] = 100;
            XfineGrid["Be"] = 100;
            XfineGrid["B"] = 100;
            XfineGrid["C"] = 100;
            XfineGrid["N"] = 100;
            XfineGrid["O"] = 100;
            XfineGrid["F"] = 100;
            XfineGrid["Ne"] = 100;

        // orders for 2nd row elements taken from NWChem
            XfineGrid["Na"] = 125;
            XfineGrid["Mg"] = 125;
            XfineGrid["Al"] = 125;
            XfineGrid["Si"] = 125;
            XfineGrid["P"] = 125;
            XfineGrid["S"] = 125;
            XfineGrid["Cl"] = 125;
            XfineGrid["Ar"] = 125;
            
        // orders for 3rd row elements taken from NWChem
            XfineGrid["K"] = 160;
            XfineGrid["Ca"] = 160;
            XfineGrid["Sc"] = 160;
            XfineGrid["Ti"] = 160;
            XfineGrid["V"] = 160;
            XfineGrid["Cr"] = 160;
            XfineGrid["Mn"] = 160;
            XfineGrid["Fe"] = 160;
            XfineGrid["Co"] = 160;
            XfineGrid["Ni"] = 160;
            XfineGrid["Cu"] = 160;
            XfineGrid["Zn"] = 160;
            XfineGrid["Ga"] = 160;
            XfineGrid["Ge"] = 160;
            XfineGrid["As"] = 160;
            XfineGrid["Se"] = 160;
            XfineGrid["Br"] = 160;
            XfineGrid["Kr"] = 160;
            
            // 4th row (selection)
            XfineGrid["Ag"] = 205 ;  
    }            
    
    inline void FillCoarseGrid(){
        
        // order for H, He (not given in NWChem, assuming same as 1st row)
            CoarseGrid["H"]  = 35;
            CoarseGrid["He"] = 35;
        
        // orders for 1st row elements taken from NWChem
            CoarseGrid["Li"] = 35;
            CoarseGrid["Be"] = 35;
            CoarseGrid["B"] = 35;
            CoarseGrid["C"] = 35;
            CoarseGrid["N"] = 35;
            CoarseGrid["O"] = 35;
            CoarseGrid["F"] = 35;
            CoarseGrid["Ne"] = 35;

        // orders for 2nd row elements taken from NWChem
            CoarseGrid["Na"] = 70;
            CoarseGrid["Mg"] = 70;
            CoarseGrid["Al"] = 70;
            CoarseGrid["Si"] = 70;
            CoarseGrid["P"] = 70;
            CoarseGrid["S"] = 70;
            CoarseGrid["Cl"] = 70;
            CoarseGrid["Ar"] = 70;
             
        // orders for 3rd row elements taken from NWChem
            CoarseGrid["K"] = 95;
            CoarseGrid["Ca"] = 95;
            CoarseGrid["Sc"] = 95;
            CoarseGrid["Ti"] = 95;
            CoarseGrid["V"] = 95;
            CoarseGrid["Cr"] = 95;
            CoarseGrid["Mn"] = 95;
            CoarseGrid["Fe"] = 95;
            CoarseGrid["Co"] = 95;
            CoarseGrid["Ni"] = 95;
            CoarseGrid["Cu"] = 95;
            CoarseGrid["Zn"] = 95;
            CoarseGrid["Ga"] = 95;
            CoarseGrid["Ge"] = 95;
            CoarseGrid["As"] = 95;
            CoarseGrid["Se"] = 95;
            CoarseGrid["Br"] = 95;
            CoarseGrid["Kr"] = 95;
                      
            // 4th row (selection)
            CoarseGrid["Ag"] = 104 ;    
    }          
    
    inline void FillXcoarseGrid(){
        
        // order for H, He (not given in NWChem, assuming same as 1st row)
            XcoarseGrid["H"]  = 21;
            XcoarseGrid["He"] = 21;
        
        // orders for 1st row elements taken from NWChem
            XcoarseGrid["Li"] = 21;
            XcoarseGrid["Be"] = 21;
            XcoarseGrid["B"] = 21;
            XcoarseGrid["C"] = 21;
            XcoarseGrid["N"] = 21;
            XcoarseGrid["O"] = 21;
            XcoarseGrid["F"] = 21;
            XcoarseGrid["Ne"] = 21;

        // orders for 2nd row elements taken from NWChem
            XcoarseGrid["Na"] = 42;
            XcoarseGrid["Mg"] = 42;
            XcoarseGrid["Al"] = 42;
            XcoarseGrid["Si"] = 42;
            XcoarseGrid["P"] = 42;
            XcoarseGrid["S"] = 42;
            XcoarseGrid["Cl"] = 42;
            XcoarseGrid["Ar"] = 42;
                     
        // orders for 3rd row elements taken from NWChem
            XcoarseGrid["K"] = 75;
            XcoarseGrid["Ca"] = 75;
            XcoarseGrid["Sc"] = 75;
            XcoarseGrid["Ti"] = 75;
            XcoarseGrid["V"] = 75;
            XcoarseGrid["Cr"] = 75;
            XcoarseGrid["Mn"] = 75;
            XcoarseGrid["Fe"] = 75;
            XcoarseGrid["Co"] = 75;
            XcoarseGrid["Ni"] = 75;
            XcoarseGrid["Cu"] = 75;
            XcoarseGrid["Zn"] = 75;
            XcoarseGrid["Ga"] = 75;
            XcoarseGrid["Ge"] = 75;
            XcoarseGrid["As"] = 75;
            XcoarseGrid["Se"] = 75;
            XcoarseGrid["Br"] = 75;
            XcoarseGrid["Kr"] = 75;
                   
            // 4th row (selection)
            XcoarseGrid["Ag"] = 84 ;  
    }            

        };

    }}
#endif	// VOTCA_XTP_EULER_MACLAURIN_H

