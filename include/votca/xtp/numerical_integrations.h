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

#ifndef __XTP_NUMERICAL_INTEGRATION__H
#define	__XTP_NUMERICAL_INTEGRATION__H

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/basisset.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/grid_containers.h>



namespace votca { namespace xtp {

    namespace ub = boost::numeric::ublas;
    
    

        class NumericalIntegration {
        public: 
            
            NumericalIntegration():density_set(false) {};

            void GridSetup(std::string type, BasisSet* bs , std::vector<QMAtom* > _atoms  );
            void FindsignificantAtoms(AOBasis* basis);
            //used for test purposes
            double StupidIntegrate( std::vector<double>& _data );
            
            std::vector<vec const *> getGridpoints();
            
            ub::matrix<double> numAOoverlap ( AOBasis* basis  );
            
            double IntegrateDensity_Atomblock(const ub::matrix<double>& _density_matrix, AOBasis* basis);
            double IntegratePotential(const vec& rvector);
            
            double IntegrateField(const std::vector<double>& externalfield);
            
            double getExactExchange(const std::string _functional);
            
            ub::matrix<double> IntegrateVXC_Atomblock (const ub::matrix<double>& _density_matrix, AOBasis* basis,const std::string _functional);
            ub::matrix<double> IntegrateExternalPotential_Atomblock(AOBasis* basis,std::vector<double> Potentialvalues);
             
            
            // this gives int (e_xc-V_xc)*rho d3r
            double getTotEcontribution(){return EXC;}
            //ub::matrix<double> StupidIntegrateVXC ( ub::matrix<double>& _density_matrix, AOBasis* basis  );
            
        private:
            
           
            
            std::vector<double> SSWpartition( int ngrid, int igrid, int ncenters ,  std::vector< std::vector<double> >& rq, double ass );
            std::vector<double> Rij;
            ub::matrix<double> Rij_mat;
            int _totalgridsize;
            double erf1c(double x);
            double erfcc(double x);
            std::vector< std::vector< GridContainers::integration_grid > > _grid;
            double EXC;
            bool density_set;
            vector< vector< vector<int> > > _significant_atoms;
            vector < int > _startIdx;
            vector < int > _blocksize;
            typedef vector< AOShell* >::iterator AOShellIterator;
            vector< vector< AOShellIterator > > _atomshells;
            vector< AOShellIterator > _singleatom;

        };

    }}
#endif	/* NUMERICAL_INTEGRATION_H */
