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
#include <boost/math/constants/constants.hpp>
#include <votca/ctp/radial_euler_maclaurin_rule.h>
#include <votca/ctp/sphere_lebedev_rule.h>
using namespace std;


namespace votca {
    namespace ctp {

        void NumericalIntegration::GridSetup(string type, BasisSet* bs, vector<QMAtom*> _atoms) {

            const double pi = boost::math::constants::pi<double>();
            // get GridContainer
            GridContainers _grids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, "medium", _grids); // this checks out 1:1 with NWChem results! AWESOME

            cout << "Radial grid summary " << endl;
            map<string, GridContainers::radial_grid>::iterator it;
            for (it = _grids._radial_grids.begin(); it != _grids._radial_grids.end(); ++it) {
                cout << " Element " << it->first << " Number of points " << it->second.radius.size() << endl;
            }

            // get angular grid per element
            LebedevGrid _sphericalgrid;
            cout << "Spherical grid summary " << endl;
            for (it = _grids._radial_grids.begin(); it != _grids._radial_grids.end(); ++it) {
                _sphericalgrid.getSphericalGrid(_atoms, "medium", _grids);
                cout << " Element " << it->first << " Number of points " << _grids._spherical_grids.at(it->first).weight.size() << endl;

            }


            // combine the element-based information with the geometry
            std::vector< GridContainers::integration_grid > _grid;
            vector< QMAtom* > ::iterator ait;
            vector< QMAtom* > ::iterator bit;
            // int i_atom = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get center coordinates
                double x_c = (*ait)->x * 1.8897259886;
                double y_c = (*ait)->y * 1.8897259886;
                double z_c = (*ait)->z * 1.8897259886;
                string name = (*ait)->type;

                // get radial grid information for this atom type
                GridContainers::radial_grid _radial_grid = _grids._radial_grids.at(name);

                // get spherical grid information for this atom type
                GridContainers::spherical_grid _spherical_grid = _grids._spherical_grids.at(name);


                // for each (theta,phi)
                for (int _i_sph = 0; _i_sph < _spherical_grid.phi.size(); _i_sph++) {

                    double p = _spherical_grid.phi[_i_sph] * pi / 180.0; // back to rad
                    double t = _spherical_grid.theta[_i_sph] * pi / 180.0; // back to rad
                    double ws = _spherical_grid.weight[_i_sph];

                    double x_s = sin(p) * cos(t);
                    double y_s = sin(p) * sin(t);
                    double z_s = cos(p);

                    // for each radial value
                    for (int _i_rad = 0; _i_rad < _radial_grid.radius.size(); _i_rad++) {
                        double r = _radial_grid.radius[_i_rad];
                        GridContainers::integration_grid _gridpoint;
                        _gridpoint.grid_x = x_c + r * x_s;
                        _gridpoint.grid_y = y_c + r * y_s;
                        _gridpoint.grid_z = z_c + r * z_s;

                        _gridpoint.grid_weight = ws * _radial_grid.weight[_i_rad];

                        //_gridpoint.grid_atom = i_atom;
                        
                        _grid.push_back(_gridpoint);


                    } // radial gridpoints
                } // spherical gridpoint
                // i_atom++;
            } // atoms


            cout << " Total size of integration grid: " << _grid.size() << endl;

            // partition function magic

            // for each center
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get center coordinates
                double x_c = (*ait)->x * 1.8897259886;
                double y_c = (*ait)->y * 1.8897259886;
                double z_c = (*ait)->z * 1.8897259886;

                std::vector<double> temp;
                // for each gridpoint
                for (std::vector<GridContainers::integration_grid >::iterator git = _grid.begin(); git != _grid.end(); ++git) {

                    double x = (*git).grid_x - x_c;
                    double y = (*git).grid_y - y_c;
                    double z = (*git).grid_z - z_c;

                    temp.push_back(sqrt(x * x + y * y + z * z));

                } // gridpoint
                rq.push_back(temp); // rq[center][gridpoint]

            } // centers

            // get all inter-center distances

            int ij = 0;
            Rij.push_back(0.0); // 1st center "self-distance"
            for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
                // get center coordinates
                double x_a = (*ait)->x * 1.8897259886;
                double y_a = (*ait)->y * 1.8897259886;
                double z_a = (*ait)->z * 1.8897259886;

                for (bit = _atoms.begin(); bit != ait; ++bit) {
                    ij++;
                    // get center coordinates
                    double x_b = (*bit)->x * 1.8897259886;
                    double y_b = (*bit)->y * 1.8897259886;
                    double z_b = (*bit)->z * 1.8897259886;

                    Rij.push_back(1.0 / sqrt((x_a - x_b)*(x_a - x_b) + (y_a - y_b)*(y_a - y_b) + (z_a - z_b)*(z_a - z_b)));


                } // atoms
                Rij.push_back(0.0); // self-distance again

            } // atoms

            // NOW the real thing... 

            // for each center again
            int i_center =0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                // find nearest neighbor
                double distNN = 1e10;

                vector< QMAtom* > ::iterator NNit;
                int i_NN;

                // get center coordinates
                double x_a = (*ait)->x * 1.8897259886;
                double y_a = (*ait)->y * 1.8897259886;
                double z_a = (*ait)->z * 1.8897259886;


                // now check all other center
                int i_b =0;
                for (bit = _atoms.begin(); bit != _atoms.end(); ++bit) {

                    if (bit != ait) {
                        // get center coordinates
                        double x_b = (*bit)->x * 1.8897259886;
                        double y_b = (*bit)->y * 1.8897259886;
                        double z_b = (*bit)->z * 1.8897259886;

                        double distSQ = (x_a - x_b)*(x_a - x_b) + (y_a - y_b)*(y_a - y_b) + (z_a - z_b)*(z_a - z_b);

                        // update NN distance and iterator
                        if ( distSQ < distNN ) {
                            distNN = distSQ;
                            NNit = bit;
                            i_NN = i_b;
                        }

                    } // if ( ait != bit) 
                    i_b++;
                }// bit centers


              
                double eps = 0.002;
                double ass = 0.725;
                double radwgh = (1.0 - ass ) * sqrt(distNN) * 0.5;
                // reduce "checklist" forward (-> get left start iterator)
                               int _left = 0;
                std::vector<GridContainers::integration_grid >::iterator leftit;
                //for (std::vector<GridContainers::integration_grid >::iterator git = _grid.begin(); git != _grid.end(); ++git) {
                for ( int i_grid  = 0 ; i_grid < _grid.size(); i_grid++) {
                    if ( rq[i_center][i_grid] > (radwgh + eps)  ) {
                        //leftit = git;
                        _left = i_grid;
                        break; // out of the for-loop
                    }
                    i_grid++;
                }
                
                // update NN distance
                distNN = (ass-eps) * sqrt(distNN) ;
                // reduce "checklist" backward (-> get right start iterator)
                //std::vector<GridContainers::integration_grid >::iterator rightit;
                //i_grid = _grid.size()-1;
                //for (std::vector<GridContainers::integration_grid >::iterator git = _grid.rbegin(); git != leftit; ++git) {
                int _right;
                for (int i_grid = _grid.size()-1; i_grid >= _left; i_grid--) {
                    if (  (rq[i_center][i_grid] - rq[i_center][i_NN] ) > distNN   ) {
                        // set weight to zero
                        _grid[i_grid].grid_weight = 0.0;
                    } else {
                        _right = i_grid;
                        break;
                    }
                } // right iterator
                
                // now for the remaining grid points
                for ( int i_grid = _left; i_grid <= _right ; i_grid++){
                    
                    
                    // call some shit called grid_ssw0 in NWChem
                    std::vector<double> _p = SSWpartition( _grid.size(), i_grid, _atoms.size(), ass );
                    // check weight sum
                    
                    
                    
                    
                    
                } // grid points
                
                // compress grid  == remove gridpoints with weight below tolerance
                
                
             
                











                
                i_center++;
            } // ait centers


        }

        
        
        std::vector<double> NumericalIntegration::SSWpartition(int ngrid, int igrid, int ncenters, double ass){
            
            // initialize partition vector to 1.0
            std::vector<double> p(ncenters,1.0);
            
            const double tol_scr = 1e-10;
            const double leps    = 1e-6; 
            // go through centers
            for ( int i = 1; i < ncenters; i++ ){
                
                int ij = i*(i+1)/2 -1; // indexing magic
                double rag = rq[i][igrid] ;
                
                // through all other centers (one-directional)
                for (int j = 0; j < i; j++ ){
                    
                    ij++;
                    if ( ( std::abs(p[i]) > tol_scr  ) || ( std::abs(p[j]) > tol_scr  ) ){
                        
                        double mu = ( rag - rq[j][igrid] )*Rij[ij]; 
                        if ( mu > ass ) {
                            p[i] = 0.0;
                        } else if ( p[j] < -ass ) {
                            p[j] = 0.0;
                        } else {
                            
                            double sk;
                            if (std::abs(mu) < leps ) {
                                sk = -1.88603178008*mu + 0.5;
                            } else {
                                //sk = erf1c(mu); // FAILS
                            }
                            if ( mu > 0.0 ) sk = 1.9 - sk;
                            p[j] = p[j] * sk;
                            p[i] = p[i] * (1.0-sk);
                            
                            
                        }
                        
                    }
                    
                    
                }
                
                
                
                
            }
            
            
            
            
            
            
        }






    }
}