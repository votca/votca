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

#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/forces.h>
#include <votca/xtp/bfgs-trm.h>

namespace votca {
    namespace xtp {
      using std::cout;
      using std::endl;
      using namespace tools;

        void BFGSTRM::Initialize(Property &options) {

            // default convergence parameters from ORCA
            _convergence = options.ifExistsReturnElseReturnDefault<double>(".convergence.energy", 1.e-6); // Hartree
            _RMSForce_convergence = options.ifExistsReturnElseReturnDefault<double>(".convergence.RMSForce", 3.e-5); // Hartree/Bohr
            _MaxForce_convergence = options.ifExistsReturnElseReturnDefault<double>(".convergence.MaxForce", 1.e-4); // Hartree/Bohr
            _RMSStep_convergence = options.ifExistsReturnElseReturnDefault<double>(".convergence.RMSStep", 6.e-4); // Bohr
            _MaxStep_convergence = options.ifExistsReturnElseReturnDefault<double>(".convergence.MaxStep", 1.e-3); // Bohr

            // initial trust radius
            _trust_radius = options.ifExistsReturnElseReturnDefault<double>(".trust", 0.01); // Angstrom

            // maximum number of iterations
            _max_iteration = options.ifExistsReturnElseReturnDefault<unsigned>(".maxiter", 50);

            // restart from saved history
            _restart_opt = options.ifExistsReturnElseReturnDefault<bool>(".restart", false);
            if (_restart_opt && !boost::filesystem::exists("optimization.restart")) {
                throw runtime_error(string("\n Restart requested but optimization.restart file not found!"));
            };


            // fix units from Ang to Bohr
            _trust_radius = _trust_radius * tools::conv::ang2bohr; // initial trust radius in a.u.
            _trust_radius_max = 0.1;
            _spintype = _force_engine.GetSpinType();
            _opt_state = _force_engine.GetOptState();

            _natoms =_orbitals.QMAtoms().size();
            _force =Eigen::MatrixX3d::Zero(_natoms,3);
            _force_old = Eigen::MatrixX3d::Zero(_natoms,3);
            _xyz_shift = Eigen::MatrixX3d::Zero(_natoms,3);
            _current_xyz = Eigen::MatrixX3d::Zero(_natoms,3);
            _old_xyz = Eigen::MatrixX3d::Zero(_natoms,3);
            _trial_xyz = Eigen::MatrixX3d::Zero(_natoms,3);
            _hessian = Eigen::MatrixXd::Zero(3 * _natoms, 3 * _natoms);

            // construct vectors (because?)
            _dim = 3 * _natoms;
            _previous_pos = Eigen::VectorXd::Zero(_dim);
            _current_pos = Eigen::VectorXd::Zero(_dim);
            _previous_gradient = Eigen::VectorXd::Zero(_dim);
            _current_gradient = Eigen::VectorXd::Zero(_dim);

            // Initial coordinates
            OrbitalsToMatrix();

            return;

        }

        void BFGSTRM::Optimize() {

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of total energy: %1$8.6f Hartree ") % _convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of RMS Force:    %1$8.6f Hartree/Bohr ") % _RMSForce_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of Max Force:    %1$8.6f Hartree/Bohr ") % _MaxForce_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of RMS Step:     %1$8.6f Bohr ") % _RMSStep_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of Max Step:     %1$8.6f Bohr ") % _MaxStep_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM initial trust radius:        %1$8.6f Bohr") % _trust_radius).str() << flush;
            bool _converged = false;

            // in 1st iteration, get the energy of the initial configuration
            if (_iteration == 0) {
                _last_energy = GetEnergy(); // in Hartree
                WriteTrajectory();
            }

            // now change geometry until convergence
            while (!_converged) {
                _iteration++;
                _update_hessian = true; // always update Hessian after each accepted step

                // calculate the forces for the present configuration
                _force_engine.Calculate(_last_energy);
                _force = _force_engine.GetForces();

                _step_accepted = false;
                while (!_step_accepted) {

                    // determine new trial configuration according to BFGS
                    BFGSStep();

                    // update coordinates in segment
                    UpdateCoordinatesOrbitals();

                    // for the updated geometry, get new reference energy
                    _new_energy = GetEnergy();
                    _energy_delta = _new_energy - _last_energy;

                    // check the energy at the trial coordinates, accept/reject step, and adjust trust radius
                    AcceptReject();


                } // checking step to be trusted

                // after the step is accepted, we can shift the stored data
                _last_energy = _new_energy;
                _old_xyz = _current_xyz;
                _current_xyz = _trial_xyz;
                _force_old = _force;

                // Write to trajectory
                WriteTrajectory();


                // Check convergence criteria
                _converged = GeometryConverged();

                // Report summary of the iteration
                Report();

                if (_iteration == _max_iteration) {
                    CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: not converged after %2$d iterations ") % _iteration % _max_iteration).str() << flush;
                    break;
                }

            } // convergence

            return;

        }

        /* Accept/reject the new geometry and adjust trust radius, if required */
        void BFGSTRM::AcceptReject() {

            if (_energy_delta > 0.0) {
                // total energy has unexpectedly increased, half the trust radius
                _trust_radius = 0.5 * _trust_radius;
                // Hessian should not be updated for reduced trust radius
                _update_hessian = false;
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: step rejected ") % _iteration).str() << flush;
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.6f Bohr") % _iteration % _trust_radius).str() << flush;
                if (_trust_radius < 1e-5) throw runtime_error("Trust radius vanishing! Stopping optimization! ");
                // repeat BGFSStep with new trust radius
            } else {
                // total energy has decreased, we accept the step but might update the trust radius
                _step_accepted = true;
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: step accepted ") % _iteration).str() << flush;
                // adjust trust radius, if required
                double _tr_check = _energy_delta / _delta_energy_estimate;
                if (_tr_check > 0.75 && 1.25 * sqrt(_norm_delta_pos) > _trust_radius) {
                    _trust_radius = 2.0 * _trust_radius;
                    // if ( _trust_radius > _trust_radius_max) _trust_radius = _trust_radius_max;
                    CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.6f Bohr") % _iteration % _trust_radius).str() << flush;
                } else if (_tr_check < 0.25) {
                    _trust_radius = 0.25 * sqrt(_norm_delta_pos);
                    CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.6f Bohr") % _iteration % _trust_radius).str() << flush;
                }

            }

            return;
        }

        /* Get the energy for the current configuration */
        double BFGSTRM::GetEnergy() {


            _gwbse_engine.setRedirectLogger(true);
            string _logger_file = "gwbse_iteration_" + (boost::format("%1%") % _iteration).str() + ".log";
            _gwbse_engine.setLoggerFile(_logger_file);
            _gwbse_engine.ExcitationEnergies(_qmpackage, _orbitals);
            double _energy = _orbitals.getTotalEnergy(_spintype, _opt_state); // in Hartree
            _gwbse_engine.setRedirectLogger(false);

            return _energy;


        }

        /* Report results of accepted step*/
        void BFGSTRM::Report() {

            // accepted step
            CTP_LOG(ctp::logINFO, *_pLog) << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format(" =========== OPTIMIZATION SUMMARY ================================= ")).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << " At iteration  " << _iteration << flush;
            _force_engine.Report();

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   ---- POSITIONS (Angstrom)   ")).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Atom\t x\t  y\t  z ")).str() << flush;
            for (unsigned _i = 0; _i < _natoms; _i++) {
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format(" %1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f")
                        % _i % (_current_xyz(_i, 0) * votca::tools::conv::bohr2ang) % (_current_xyz(_i, 1) * votca::tools::conv::bohr2ang) % (_current_xyz(_i, 2) * votca::tools::conv::bohr2ang)).str() << flush;
            }

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Total energy:     %1$12.8f Hartree ") % _new_energy).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   quadratic energy: %1$12.8f Hartree ") % _delta_energy_estimate).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   energy change:    %1$12.8f Hartree      %2$s (%3%)") % _energy_delta % Converged(_energy_converged) % _convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   RMS force:        %1$12.8f Hartree/Bohr %2$s (%3%)") % _RMSForce % Converged(_RMSForce_converged) % _RMSForce_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Max force:        %1$12.8f Hartree/Bohr %2$s (%3%)") % _MaxForce % Converged(_MaxForce_converged) % _MaxForce_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   RMS step:         %1$12.8f Bohr         %2$s (%3%)") % _RMSStep % Converged(_RMSStep_converged) % _RMSStep_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Max step:         %1$12.8f Bohr         %2$s (%3%)") % _MaxStep % Converged(_MaxStep_converged) % _MaxStep_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Trust radius:     %1$12.8f Bohr     ") % _trust_radius).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << flush;
            return;
        }

        /* Convergence output */
        string BFGSTRM::Converged(bool converged) {

            if (converged) {
                return "converged    ";
            } else {
                return "not converged";
            }

        }

        /* Check convergence */
        bool BFGSTRM::GeometryConverged() {


            // checking convergence
            _energy_converged = false;
            _RMSForce_converged = false;
            _MaxForce_converged = false;
            _RMSStep_converged = false;
            _MaxStep_converged = false;

            _xyz_shift = _current_xyz - _old_xyz;

            _RMSForce = _force.cwiseAbs2().sum()/(_force.rows()*_force.cols());
            _MaxForce = _force.cwiseAbs().maxCoeff();
            _RMSStep = _xyz_shift.cwiseAbs2().sum()/(_xyz_shift.rows()*_xyz_shift.cols());
            _MaxStep = _xyz_shift.cwiseAbs().maxCoeff();

            if (std::abs(_energy_delta) < _convergence) _energy_converged = true;
            if (std::abs(_RMSForce) < _RMSForce_convergence) _RMSForce_converged = true;
            if (std::abs(_MaxForce) < _MaxForce_convergence) _MaxForce_converged = true;
            if (std::abs(_RMSStep) < _RMSStep_convergence) _RMSStep_converged = true;
            if (std::abs(_MaxStep) < _MaxStep_convergence) _MaxStep_converged = true;

            if (_energy_converged && _RMSForce_converged && _MaxForce_converged && _RMSStep_converged && _MaxStep_converged) {
                return true;
            } else {
                return false;
            }


        }

        void BFGSTRM::UpdateCoordinatesOrbitals() {
            std::vector<QMAtom*>& atoms= _orbitals.QMAtoms();
            for ( int i = 0; i< atoms.size(); i++) {
                tools::vec pos_displaced(_trial_xyz(i, 0), _trial_xyz(i, 1) , _trial_xyz(i, 2));
                atoms[i]->setPos(pos_displaced);
            }
            return;

        }

        /* Determine new trial coordinates according to BFGS */
        void BFGSTRM::BFGSStep() {

            // Rewrite2Vectors (let's rethink that later)
            Rewrite2Vectors();

            // Update Hessian
            if (_update_hessian) UpdateHessian();

            // Get displacement of coordinates
            PredictDisplacement();

            // TRM -> trust radius check and regularization, if needed
            if (OutsideTrustRegion(_norm_delta_pos)) RegularizeStep();

            // expected energy change on quadratic surface
            QuadraticEnergy();

            // new trial coordinated are written to _trial_xyz
            Rewrite2Matrices();

            return;
        }

        /* Re-Store the matrix data into vectors (let's rethink this later) */
        void BFGSTRM::Rewrite2Vectors() {
            for (unsigned _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (unsigned _i_cart = 0; _i_cart < 3; _i_cart++) {

                    int _idx = 3 * _i_atom + _i_cart;
                    _previous_pos(_idx) = _old_xyz(_i_atom, _i_cart);
                    _current_pos(_idx) = _current_xyz(_i_atom, _i_cart);
                    _previous_gradient(_idx) = -_force_old(_i_atom, _i_cart);
                    _current_gradient(_idx) = -_force(_i_atom, _i_cart);

                }
            }

            return;

        }

        /* Update the Hessian */
        void BFGSTRM::UpdateHessian() {

            // delta is new - old
            _delta_pos = _current_pos - _previous_pos;
            _norm_delta_pos = _delta_pos.squaredNorm();

            // we have no Hessian in the first iteration => start with something
            if (_iteration == 1) {
             _hessian=Eigen::MatrixXd::Identity(_dim,_dim);
            }else {
                /* for later iteration, we can make use of an iterative refinement of
                 * the initial Hessian based on the gradient (force) history
                 */

                Eigen::VectorXd _delta_gradient = _current_gradient - _previous_gradient;

                // second term in BFGS update (needs current Hessian)
               _hessian -= _hessian*_delta_pos*_delta_pos.transpose()*_hessian.transpose() / (_delta_pos.transpose()*_hessian*_delta_pos).value();

                // first term in BFGS update
                _hessian += (_delta_gradient* _delta_gradient.transpose()) / (_delta_gradient.transpose()*_delta_pos);

                // symmetrize Hessian (since d2E/dxidxj should be symmetric)
               _hessian=0.5*(_hessian+_hessian.transpose());
            } // update Hessian

            return;
        }

        /* Predict displacement of atom coordinates */
        void BFGSTRM::PredictDisplacement() {

            _delta_pos=_hessian.colPivHouseholderQr().solve(-_current_gradient);
            // new displacements for the atoms
            _norm_delta_pos = _delta_pos.squaredNorm();

            return;
        }

        /* Check if predicted displacement leaves trust region */
        bool BFGSTRM::OutsideTrustRegion(const double& _step) {

            double _trust_radius_squared = _trust_radius * _trust_radius;
            if (_step > _trust_radius_squared) {
                return true;
            } else {
                return false;
            }

        }

        /* Regularize step in case of prediction outside of Trust Region */
        void BFGSTRM::RegularizeStep() {

            double _max_step_squared = 0.0;

            // get eigenvalues and eigenvectors of Hessian

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_hessian);



            // start value for lambda  a bit lower than lowest eigenvalue of Hessian
            double _lambda;
            if (es.eigenvalues()(0) > 0.0) {
                _lambda = -0.05 * std::abs(es.eigenvalues()(0));
            } else {
                _lambda = 1.05 * es.eigenvalues()(0);
            }

            // for constrained step, we expect
            _max_step_squared = _norm_delta_pos;
            while (OutsideTrustRegion(_max_step_squared)) {
                _max_step_squared = 0.0;
                _lambda -= 0.05 * std::abs(es.eigenvalues()(0));
                for (unsigned _i = 0; _i < _dim; _i++) {
                    double _temp = es.eigenvectors().col(_i).transpose()* _current_gradient;
                    _max_step_squared += _temp * _temp / (es.eigenvalues()(_i) - _lambda) / (es.eigenvalues()(_i) - _lambda);
                }
            }

            _delta_pos = Eigen::VectorXd::Zero(_dim);
            for (unsigned _i = 0; _i < _dim; _i++) {
                _delta_pos -= es.eigenvectors().col(_i) * (es.eigenvectors().col(_i).transpose()*_current_gradient) / (es.eigenvalues()(_i) - _lambda);
            }

            _norm_delta_pos = _delta_pos.squaredNorm(); //_max_step_squared;

            return;
        }

        /* Estimate energy change based on quadratic approximation */
        void BFGSTRM::QuadraticEnergy() {
            _delta_energy_estimate = (_current_gradient.transpose()* _delta_pos).value() + (0.5 * _delta_pos.transpose()*_hessian*_delta_pos).value();
            return;
        }

        /* Rewrite the vector data back to matrices (to rethink) */
        void BFGSTRM::Rewrite2Matrices() {
            Eigen::VectorXd _new_pos = _current_pos + _delta_pos;
            Eigen::VectorXd _total_shift=Eigen::VectorXd::Zero(3);
            for (unsigned _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (unsigned _i_cart = 0; _i_cart < 3; _i_cart++) {
                    unsigned _idx = 3 * _i_atom + _i_cart;
                    _trial_xyz(_i_atom, _i_cart) = _new_pos(_idx);
                    _total_shift(_i_cart) += _delta_pos(_idx);
                }
            }
            return;
        }

        void BFGSTRM::OrbitalsToMatrix() {

           std::vector<QMAtom*> atoms = _orbitals.QMAtoms();
           for (unsigned i=0;i<atoms.size();i++) {
                // put trial coordinates (_current_xyz is in Bohr, segments in nm)
                _current_xyz(i, 0) = atoms[i]->getPos().getX();
                _current_xyz(i, 1) = atoms[i]->getPos().getY();
                _current_xyz(i, 2) = atoms[i]->getPos().getZ();
            }

            return;



        }

        /* Write accepted geometry to xyz trajectory file */
        void BFGSTRM::WriteTrajectory() {

           std::vector< QMAtom* > _atoms = _orbitals.QMAtoms();
            // write logger to log file
            ofstream ofs;
            string _trajectory_file = "optimization.trj";
            if (_iteration == 0) {
                ofs.open(_trajectory_file.c_str(), ofstream::out);
            } else {
                ofs.open(_trajectory_file.c_str(), ofstream::app);
            }

            if (!ofs.is_open()) {
                throw runtime_error("Bad file handle: " + _trajectory_file);
            }

            // write coordinates as xyz file
            ofs << _atoms.size() << endl;
            ofs << "iteration " << _iteration << " energy " << _last_energy << " Hartree" << endl;

            for (const QMAtom* atom:_atoms) {
                tools::vec pos=atom->getPos()*tools::conv::bohr2ang;
                // put trial coordinates (_current_xyz is in Bohr, segments in nm)
                ofs << atom->getType() << " " << pos.getX() << " " << pos.getY() << " " << pos.getZ() << endl;
            }
            ofs.close();
            return;


        }

    }
}
