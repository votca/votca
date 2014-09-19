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

#ifndef __CTP_EULER_MACLAURIN__H
#define	__CTP_EULER_MACLAURIN__H


#include <votca/tools/property.h>
#include <votca/ctp/basisset.h>
#include <votca/ctp/qmatom.h>
#include <votca/ctp/grid_containers.h>
using namespace std;


namespace votca { namespace ctp {


        class EulerMaclaurinGrid {
        public: 
            
            EulerMaclaurinGrid() { FillGrids(); };
            
            void getRadialGrid( BasisSet* bs , vector<QMAtom* > _atoms , string type, GridContainers& _grids );



        private:
            
            struct min_exp {
                double alpha;
                int l;
                double range;
            };
            
            struct grid_element {
                std::vector<double> gridpoint;
                std::vector<double> weight;
            };
       
            map<string,min_exp> _element_ranges;
            map<string,grid_element> _element_grids;
            
            int getGrid(string element, string type);
            
            double DetermineCutoff( double alpha, int l, double eps );
            double getNeglected( double alpha, int l, double cutoff);
            double RadialIntegral(double alpha, int l, double cutoff);
            
            void getRadialCutoffs( vector<QMAtom* > _atoms ,  BasisSet* bs ,string gridtype );
            void setGrid(int numberofpoints, double cutoff, std::vector<double>& point, std::vector<double>& weight );
            
            std::map<std::string, int>    MediumGrid;
            std::map<std::string, double> Accuracy;
            
            inline void FillGrids(){
        
                FillAccuracy();
                FillMediumGrid();
        
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
            
            
    }            
            
            

        };

    }}
#endif	/* EULER_MACLAURIN_H */