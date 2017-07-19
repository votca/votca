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

#ifndef __XTP_NUMERICAL_INTEGRATION__H
#define	__XTP_NUMERICAL_INTEGRATION__H

#ifdef LIBXC
#include <xc.h>
#undef LOG
#endif

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/basisset.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/grid_containers.h>
#include <votca/xtp/vxc_functionals.h>
#include <votca/xtp/exchange_correlation.h>
#include <votca/xtp/gridbox.h>
#include <votca/ctp/qmatom.h>


namespace votca { namespace xtp {

    namespace ub = boost::numeric::ublas;
    
    

        class NumericalIntegration {
        public: 
            
            NumericalIntegration():density_set(false),setXC(false) {};
            
            
            ~NumericalIntegration(){};
            
            void GridSetup(std::string type, BasisSet* bs , std::vector<ctp::QMAtom* > _atoms,AOBasis* basis  );
            
         
            double getExactExchange(const std::string _functional);
            std::vector<const vec*> getGridpoints();
            
            unsigned getGridSize() const{return _totalgridsize;}
            unsigned getBoxesSize() const{return _grid_boxes.size();}
            
            void setXCfunctional(const string _functional);
            
            double IntegrateDensity(const ub::matrix<double>& _density_matrix);
            double IntegratePotential(const vec& rvector);
            double IntegrateField(const std::vector<double>& externalfield);
            ub::matrix<double> IntegrateExternalPotential(const std::vector<double>& Potentialvalues);
            
            ub::vector<double> IntegrateGyrationTensor(const ub::matrix<double>& _density_matrix);
            
           
           
            ub::matrix<double> IntegrateVXC (const ub::matrix<double>& _density_matrix);
            
           
            
            // this gives int (e_xc-V_xc)*rho d3r
            double getTotEcontribution(){return EXC;}
          
            
        private:
            
            
           void FindSignificantShells();
            
           void EvaluateXC(const double rho,const ub::matrix<double>& grad_rho,double& f_xc, double& df_drho, double& df_dsigma);
          
           
           
           double erf1c(double x);
           double erfcc(double x);
           std::vector<double> SSWpartition(int igrid, int ncenters ,  std::vector< std::vector<double> >& rq );
           void SortGridpointsintoBlocks(std::vector< std::vector< GridContainers::integration_grid > >& grid);
            
            std::vector<double> Rij;
            AOBasis* _basis;

            double  _totalgridsize;
            
            std::vector< GridBox > _grid_boxes;
            
            
            ExchangeCorrelation _xc;
            bool _use_votca;
            int xfunc_id;
            
            
            
            
            double EXC;
            bool density_set;
            bool setXC;
            
            
            #ifdef LIBXC
            bool _use_separate;
            int cfunc_id;
            xc_func_type xfunc; // handle for exchange functional
            xc_func_type cfunc; // handle for correlation functional
            #endif
            
        };

    }}
#endif	/* NUMERICAL_INTEGRATION_H */
