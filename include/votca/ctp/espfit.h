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

#ifndef __CTP_ESPFIT__H
#define	__CTP_ESPFIT__H


#include <votca/ctp/elements.h>
#include <votca/ctp/grid.h>
#include <votca/ctp/aobasis.h>


/**
* \brief Takes a list of atoms, and the corresponding density matrix and puts out a table of partial charges
*
* 
* 
*/
using namespace std;
using namespace votca::tools;


namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
class Espfit{
public:
    
    Espfit(){};
   ~Espfit(){};
    
   void setLog(Logger *log) { _log = log; }
    
    void EvaluateAPECharges(Grid& _targetgrid, Grid& _chargepositions);
  
    void FitAPECharges(Grid& _targetgrid_fg, Grid& _targetgrid_bg, Grid& _chargepositions, double& netcharge);
    // on grid very fast
    void Fit2Density(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis,BasisSet &bs,string gridsize, bool _do_transition, double _netcharge=0.0);
    // not so fast
    void Fit2Density_analytic(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis, bool _do_transition, double _netcharge=0.0);
private:
    
     Logger *_log;
     Elements _elements; 
     string _method;
     
 
    ub::vector<double> EvalNuclearPotential( vector< QMAtom* >& _atoms, Grid _grid );
   
     // Fits partial charges to Potential on a grid, constrains net charge
    std::vector<double> FitPartialCharges( std::vector< ub::vector<double> >& _fitcenters, Grid& _grid, ub::vector<double>& _potential, double& _netcharge );
    
};
}}

#endif	/* ESPFIT_H */
