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

#ifndef __XTP_ESPFIT__H
#define	__XTP_ESPFIT__H


#include <votca/xtp/elements.h>
#include <votca/xtp/grid.h>
#include <votca/xtp/aobasis.h>
#include <votca/ctp/apolarsite.h>

/**
* \brief Takes a list of atoms, and the corresponding density matrix and puts out a table of partial charges
*
* 
* 
*/
using namespace votca::tools;


namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
class Espfit{
public:
    
    Espfit(ctp::Logger *log):_do_Transition(false),_do_svd(false) {_log = log;}
   ~Espfit(){};
   
   void setUseSVD(bool do_svd,double conditionnumber){_do_svd=do_svd;_conditionnumber=conditionnumber;}
    
    
    // on grid very fast
    void Fit2Density(std::vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis,std::string gridsize);
    // not so fast
    void Fit2Density_analytic(std::vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis);
private:
    
     ctp::Logger *_log;
     Elements _elements; 
     bool _do_Transition;
     bool _do_svd;
     double _conditionnumber;
     
     
    double getNetcharge( std::vector< QMAtom* >& _atoms, double N );
 
    ub::vector<double> EvalNuclearPotential( std::vector< QMAtom* >& _atoms, Grid _grid );
   
     // Fits partial charges to Potential on a grid, constrains net charge
    std::vector<double> FitPartialCharges( std::vector< tools::vec >& _fitcenters, Grid& _grid, ub::vector<double>& _potential, double& _netcharge );
    
};
}}

#endif	/* ESPFIT_H */
