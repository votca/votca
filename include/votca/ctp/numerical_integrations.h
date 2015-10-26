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

#ifndef __CTP_NUMERICAL_INTEGRATION__H
#define	__CTP_NUMERICAL_INTEGRATION__H

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/ctp/basisset.h>
#include <votca/ctp/aobasis.h>
#include <votca/ctp/grid_containers.h>
using namespace std;


namespace votca { namespace ctp {

    namespace ub = boost::numeric::ublas;
    
    

        class NumericalIntegration {
        public: 
            
            NumericalIntegration():density_set(false) {};

            void GridSetup(string type, BasisSet* bs , vector<QMAtom* > _atoms  );

            double StupidIntegrate( std::vector<double>& _data );
            
            void getGridpoints( ub::matrix<double>& _gridpoints );
            
            ub::matrix<double> numAOoverlap ( AOBasis* basis  );
            double IntegrateDensity(ub::matrix<double>& _density_matrix, AOBasis* basis);
            double IntegrateDensity_Atomblock(ub::matrix<double>& _density_matrix, AOBasis* basis);
            double IntegratePotential(ub::vector<double> rvector);
            ub::matrix<double> IntegrateVXC ( ub::matrix<double>& _density_matrix, AOBasis* basis  );
            ub::matrix<double> IntegrateVXC_block ( ub::matrix<double>& _density_matrix, AOBasis* basis   );
            ub::matrix<double> IntegrateVXC_Atomblock ( ub::matrix<double>& _density_matrix, AOBasis* basis );
            
            // this gives int (e_xc-V_xc)*rho d3r
            double& getTotEcontribution(){return EXC;}
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

        };

    }}
#endif	/* NUMERICAL_INTEGRATION_H */