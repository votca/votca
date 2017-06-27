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

#include <votca/tools/linalg.h>
#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/forces.h>

namespace votca {
    namespace xtp {

        void GeometryOptimization::Initialize(Property *options) {

            // checking if there is only one segment
            _nsegments = _segments.size();
            if (_nsegments > 1) throw runtime_error(string("\n Geometry optimization of more than 1 conjugated segment not supported. Stopping!"));
            
            _opt_state    = options->ifExistsReturnElseReturnDefault<int>(".state", 1);
            _spintype     = options->ifExistsReturnElseReturnDefault<string>(".spintype", "singlet");
            _convergence  = options->ifExistsReturnElseReturnDefault<double>(".convergence", 0.01); // units?

            // pre-check optimizer method
            std::vector<string> choices={ "BFGS-TRM" };
            _optimizer         = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(".optimizer.method", choices);
            _optimizer_options = options->get(".optimizer");

            // pre-check forces method
            choices={ "forward", "central" };
            _force_method  = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(".forces.method", choices);
            // _displacement = options->get(".displacement").as<double> ();
            _force_options = options->get(".forces");
            
            // BFGS-TRM specifics, maybe generalize
            _trust_radius = options->get(".optimizer.trust").as<double> ();
            
            // restart from saved history
            _restart_opt = options->ifExistsReturnElseReturnDefault<bool>(".restart", false);
            if (_restart_opt && !boost::filesystem::exists("optimization.restart")) {
                throw runtime_error(string("\n Restart requested but optimization.restart file not found!")); 
            };
            

            
            /* _trust_radius = _trust_radius * tools::conv::ang2bohr; // initial trust radius in a.u.
             _trust_radius_max = 0.1;
             _displacement = _displacement * tools::conv::ang2bohr; // initial trust radius in a.u.
             _pLog.setReportLevel(ctp::logINFO);
             CTP_LOG(ctp::logINFO, _log) << "Requested geometry optimization of excited state " << _spintype << " " << _opt_state << flush;
             CTP_LOG(ctp::logINFO, _log) << "... forces calculation using " << _forces << " differences " << flush;
             CTP_LOG(ctp::logINFO, _log) << "... displacement for numerical forces:           " << _displacement << " Bohr" << flush;
             CTP_LOG(ctp::logINFO, _log) << "... convergence of total energy:                 " << _convergence << " Hartree" << flush;
             CTP_LOG(ctp::logINFO, _log) << "... initial trust radius:                        " << _trust_radius << " Bohr " << flush;
             _pLog.setReportLevel(_ReportLevel);*/
            


            
            _natoms = _segments[0]->Atoms().size();
            _force = ub::zero_matrix<double>(_natoms, 3);
            _force_old = ub::zero_matrix<double>(_natoms, 3);
            _xyz_shift = ub::zero_matrix<double>(_natoms, 3);
            _current_xyz = ub::zero_matrix<double>(_natoms, 3);
            _old_xyz = ub::zero_matrix<double>(_natoms, 3);
            _trial_xyz = ub::zero_matrix<double>(_natoms, 3);
            _hessian = ub::zero_matrix<double>(3 * _natoms, 3 * _natoms);
            return;

        }

        void GeometryOptimization::Evaluate() {

            _gwbse_engine.ExcitationEnergies(_qmpackage, _segments, _orbitals);
            
            ctp::TLogLevel _ReportLevel = _pLog->getReportLevel(); // backup report level
            _pLog->setReportLevel(ctp::logERROR);
            
            
            double _last_energy = _orbitals->GetTotalEnergy(_spintype, _opt_state); // in Hartree
            
            cout << " Last Energy " << _last_energy << endl;
                        
            // get a force object
            Forces _force_engine(_gwbse_engine,_qmpackage, _segments, _orbitals);
            _force_engine.Initialize( &_force_options );
            _force_engine.SetOptState(_opt_state);
            _force_engine.SetSpinType( _spintype );


            // Actually calculate the forces
            _force_engine.Calculate( _last_energy );
            
            // Let's see the forces
            _force = _force_engine.GetForces();
            for ( int i = 0; i < _force.size1(); i++ ){
                for ( int j = 0; j< _force.size2(); j++ ){
                    
                    cout << " atom " << i << " direction " << j << " force: " << _force(i,j) << "\n" << endl;
                    
                }
                
            }
            
            _pLog->setReportLevel( _ReportLevel ); // 
            
            
            
            
            
            


        }

        void GeometryOptimization::Checkpoint(std::vector<ctp::Segment*>& _molecule) {

            // writing current coordinates and forces
            FILE *out = NULL;
            std::vector<ctp::Segment*> ::iterator sit;
            string FILENAME = "geometry_optimization.xyz";
            for (sit = _molecule.begin(); sit < _molecule.end(); ++sit) {
                if (_iteration == 0) {
                    out = fopen(FILENAME.c_str(), "w");
                } else {
                    out = fopen(FILENAME.c_str(), "a");
                }
                WriteIteration(out, (*sit));
            }

            fclose(out);

        }



        // write iteration

        void GeometryOptimization::WriteIteration(FILE *out, ctp::Segment* _segment) {

            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Atom* > _atoms = _segment->Atoms();
            int qmatoms = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                ++qmatoms;
            }

            fprintf(out, "%6d \n", qmatoms);
            fprintf(out, "at iteration %d \n", _iteration);

            int _i_atom = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                votca::tools::vec pos = (*ait)->getQMPos();
                string name = (*ait)->getElement();

                fprintf(out, "%2s %4.7f %4.7f %4.7f atom_vector %4.7f %4.7f %4.7f \n",
                        name.c_str(),
                        pos.getX()*10,
                        pos.getY()*10,
                        pos.getZ()*10,
                        _force(_i_atom, 0)*51.4220822067,
                        _force(_i_atom, 1)*51.4220822067,
                        _force(_i_atom, 2)*51.4220822067);
                _i_atom++;
                return;
            }

            return;
        }

        void GeometryOptimization::BFGSStep(int& _iteration, bool& _update_hessian, ub::matrix<double>& _force, ub::matrix<double>& _force_old, ub::matrix<double>& _current_xyz, ub::matrix<double>& _old_xyz, ub::matrix<double>& _hessian, ub::matrix<double>& _shift, ub::matrix<double>& _trial_xyz) {

            // stuff to prepare
            ub::vector<double> _total_force(3, 0.0);
            double _force_norm = 0.0;
            for (int _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (int _i_cart = 0; _i_cart < 3; _i_cart++) {
                    _total_force(_i_cart) += _force(_i_atom, _i_cart);
                    _force_norm += _force(_i_atom, _i_cart) * _force(_i_atom, _i_cart);
                }
            }

            cout << "Total force: " << _total_force(0) << " " << _total_force(1) << " " << _total_force(2) << endl;

            // remove CoM force
            //CTP_LOG(ctp::logINFO,_log) << " ----------------- Forces ------------------- " << flush;
            //CTP_LOG(ctp::logINFO,_log) << " Atom        x-comp       y-comp      z-comp  " << flush;
            for (int _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (int _i_cart = 0; _i_cart < 3; _i_cart++) {
                    // _force(_i_atom,_i_cart) -= _total_force(_i_cart)/_natoms;
                }
                //  CTP_LOG(ctp::logINFO,_log) <<  _i_atom << "   " <<  _force(_i_atom,0) << "   " <<  _force(_i_atom,1) << "   " <<  _force(_i_atom,2) << endl;
            }
            //CTP_LOG(ctp::logINFO,_log) << " -------------------------------------------- " << flush;



            // construct vectors
            int _dim = 3 * _natoms;
            ub::vector<double> _previous_pos(_dim, 0.0);
            ub::vector<double> _current_pos(_dim, 0.0);
            ub::vector<double> _previous_gradient(_dim, 0.0);
            ub::vector<double> _current_gradient(_dim, 0.0);

            for (int _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (int _i_cart = 0; _i_cart < 3; _i_cart++) {

                    int _idx = 3 * _i_atom + _i_cart;
                    _previous_pos(_idx) = _old_xyz(_i_atom, _i_cart);
                    _current_pos(_idx) = _current_xyz(_i_atom, _i_cart);
                    _previous_gradient(_idx) = -_force_old(_i_atom, _i_cart);
                    _current_gradient(_idx) = -_force(_i_atom, _i_cart);

                }
            }

            // delta is new - old
            ub::vector<double> _delta_pos = _current_pos - _previous_pos;
            _norm_delta_pos = ub::inner_prod(_delta_pos, _delta_pos);


            // update of Hessian only in first trial 
            if (_update_hessian) {

                // we have no Hessian in the first iteration => start with something
                if (_iteration == 0) {

                    for (int _i = 0; _i < _dim; _i++) {
                        _hessian(_i, _i) = 1.0; // unit matrix
                    }

                } else {
                    /* for later iteration, we can make use of an iterative refinement of 
                     * the initial Hessian based on the gradient (force) history
                     */

                    ub::vector<double> _delta_gradient = _current_gradient - _previous_gradient;

                    // second term in BFGS update (needs current Hessian)
                    ub::vector<double> _temp1 = ub::prod(_hessian, _delta_pos);
                    ub::vector<double> _temp2 = ub::prod(ub::trans(_delta_pos), _hessian);
                    _hessian -= ub::outer_prod(_temp1, _temp2) / ub::inner_prod(_delta_pos, _temp1);

                    // first term in BFGS update
                    _hessian += ub::outer_prod(_delta_gradient, _delta_gradient) / ub::inner_prod(_delta_gradient, _delta_pos);

                    // symmetrize Hessian (since d2E/dxidxj should be symmetric)
                    for (int _i = 0; _i < _dim; _i++) {
                        for (int _j = _i + 1; _j < _dim; _j++) {
                            double _sum = 0.5 * (_hessian(_i, _j) + _hessian(_j, _i));
                            _hessian(_i, _j) = _sum;
                            _hessian(_j, _i) = _sum;
                        }
                    }

                }
            } // update Hessian


            // get inverse of the Hessian
            ub::matrix<double> _hessian_inverse;
            votca::tools::linalg_invert(_hessian, _hessian_inverse);

            // new displacements for the atoms
            _delta_pos = -ub::prod(_hessian_inverse, _current_gradient);
            _norm_delta_pos = ub::inner_prod(_delta_pos, _delta_pos);

            // TRM -> trust radius check
            double _trust_radius_squared = _trust_radius * _trust_radius;
            double _max_step_squared = 0.0;

            // cout << " step sq " << _norm_delta_pos << " vs. TR sq" << _trust_radius_squared << endl;
            if (_norm_delta_pos > _trust_radius_squared) {

                // get eigenvalues and eigenvectors of Hessian
                ub::matrix<double> _eigenvectors;
                ub::vector<double> _eigenvalues;
                votca::tools::linalg_eigenvalues(_hessian, _eigenvalues, _eigenvectors);

                // CTP_LOG(ctp::logDEBUG, _log) << " BFGS-TRM: Lowest eigenvalue of Hessian is... " << _eigenvalues(0) << flush;

                // start value for lambda  a bit lower than lowest eigenvalue of Hessian
                double _lambda;
                if (_eigenvalues(0) > 0.0) {
                    _lambda = -0.05 * std::abs(_eigenvalues(0));
                } else {
                    _lambda = 1.05 * _eigenvalues(0);
                }

                // for constrained step, we expect
                _max_step_squared = _norm_delta_pos;
                while (_max_step_squared > _trust_radius_squared) {
                    _max_step_squared = 0.0;
                    _lambda -= 0.05 * std::abs(_eigenvalues(0));
                    for (int _i = 0; _i < _dim; _i++) {
                        ub::vector<double> _slice(_dim, 0.0);
                        for (int _j = 0; _j < _dim; _j++) {
                            _slice(_j) = _eigenvectors(_j, _i);
                        }

                        //cout << " forece is of dim " << _current_force.size1() << "  " << _current_force.size2() << endl;
                        double _temp = ub::inner_prod(_slice, _current_gradient);
                        //cout << " slice is of dim " << _slice.size1() << "  " << _slice.size2() << endl;            cout << " tmep is of dim " << _temp.size1() << "  " << _temp.size2() << endl;
                        // cout << " into max_step_sq " << _temp << " and  "  << ( _eigenvalues(_i) - _lambda ) << endl;
                        _max_step_squared += _temp * _temp / (_eigenvalues(_i) - _lambda) / (_eigenvalues(_i) - _lambda);
                    }
                }

                //CTP_LOG(ctp::logDEBUG,_log) << " BFGS-TRM: with lambda " << _lambda << " max step sq is " << _max_step_squared << flush;

                _delta_pos = ub::zero_vector<double>(_dim);
                for (int _i = 0; _i < _dim; _i++) {
                    ub::vector<double> _slice(_dim, 0.0);
                    for (int _j = 0; _j < _dim; _j++) {
                        _slice(_j) = _eigenvectors(_j, _i);
                    }

                    _delta_pos -= _slice * ub::inner_prod(_slice, _current_gradient) / (_eigenvalues(_i) - _lambda);
                }

                _norm_delta_pos = ub::inner_prod(_delta_pos, _delta_pos); //_max_step_squared; 
            }






            ub::vector<double> _new_pos = _current_pos + _delta_pos;
            //CTP_LOG(ctp::logDEBUG,_log) << "BFGS-TRM: step " << sqrt(_norm_delta_pos) << " vs TR " << sqrt(_trust_radius_squared) << flush  ;

            // get the energy estimate for a local quadratic surface
            ub::vector<double> _temp_ene = ub::prod(_hessian, _delta_pos);
            _delta_energy_estimate = ub::inner_prod(_current_gradient, _delta_pos) + 0.5 * ub::inner_prod(_delta_pos, _temp_ene);

            //CTP_LOG(ctp::logINFO,_log) << " BFGS-TRM: estimated energy change: " << setprecision(12) << _delta_energy_estimate << endl;

            // update atom coordinates
            ub::vector<double> _total_shift(3, 0.0);
            for (int _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (int _i_cart = 0; _i_cart < 3; _i_cart++) {
                    int _idx = 3 * _i_atom + _i_cart;
                    //_current_xyz(_i_atom,_i_cart) = _new_pos(_idx);
                    _trial_xyz(_i_atom, _i_cart) = _new_pos(_idx);
                    _total_shift(_i_cart) += _delta_pos(_idx);
                }
            }
            /* 
             // make sure there is no CoM movement
                for ( int _i_atom = 0; _i_atom < _natoms; _i_atom++){
                 for ( int _i_cart = 0; _i_cart < 3; _i_cart++ ){
                      //int _idx = 3*_i_atom + _i_cart;
                     //_current_xyz(_i_atom,_i_cart) -= _total_shift(_i_cart)/_natoms;
                     // _trial_xyz(_i_atom,_i_cart) -= _total_shift(_i_cart)/_natoms;
                 }
             }
             */
            return;
        }


    }
}
