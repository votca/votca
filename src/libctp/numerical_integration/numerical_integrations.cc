/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICEN_olE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>


#include <votca/ctp/numerical_integrations.h>
#include <votca/ctp/radial_euler_maclaurin_rule.h>
#include <votca/ctp/sphere_lebedev_rule.h>
using namespace std;


namespace votca { namespace ctp {

        void NumericalIntegration::GridSetup(string type, BasisSet* bs, vector<QMAtom*> _atoms) {

            // get GridContainer
            GridContainers _grids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, "medium", _grids); // this checks out 1:1 with NWChem results! AWESOME

            cout << "Radial grid summary " << endl;
            map<string, GridContainers::radial_grid>::iterator it;
            for ( it = _grids._radial_grids.begin() ; it != _grids._radial_grids.end() ; ++it ){
                
                cout << " Element " << it->first << " Number of points " << it->second.radius.size() <<  endl;
                
            }

            // get angular grid per element
            LebedevGrid _sphericalgrid;
                cout << "Spherical grid summary " << endl;
            for ( it = _grids._radial_grids.begin() ; it != _grids._radial_grids.end() ; ++it ){
                _sphericalgrid.getSphericalGrid( _atoms, "medium", _grids );
                
                
                  cout << " Element " << it->first << " Number of points " << _grids._spherical_grids.at(it->first).weight.size() <<  endl;
                
            }
            
            


      /*     // try getting a Lebedev Grid

           
            std::vector<double> _theta;
            std::vector<double> _phi;
            std::vector<double> _weight;
 
            _lebedevgrid.getUnitSphereGrid("C","medium", _theta, _phi, _weight);
            
            
            cout << endl;
            double _wsum =0.0;
            for ( int _i =0; _i < _theta.size() ; _i++ ){
                 // _lebedevgrid.xyz_to_tp(x[_i],y[_i],z[_i],&_theta,&_phi);
                cout << _i << ":(" <<  _theta[_i]  << "," << _phi[_i] << ") = = " << _weight[_i] << endl;
                _wsum += _weight[_i];
            }

            cout << _wsum << endl;
            exit(0);*/

            
            // combine


        }







}
}