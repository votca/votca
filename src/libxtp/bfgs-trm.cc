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

#include "votca/xtp/statefilter.h"

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
            _trust_radius = _trust_radius * tools::conv::ang2bohr; // initial trust radius in a.u.

            _max_iteration = options.ifExistsReturnElseReturnDefault<unsigned>(".maxiter", 50);

            return;
        }

        void BFGSTRM::Optimize() {

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of total energy: %1$8.6f Hartree ") % _convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of RMS Force:    %1$8.6f Hartree/Bohr ") % _RMSForce_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of Max Force:    %1$8.6f Hartree/Bohr ") % _MaxForce_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of RMS Step:     %1$8.6f Bohr ") % _RMSStep_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM convergence of Max Step:     %1$8.6f Bohr ") % _MaxStep_convergence).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM initial trust radius:        %1$8.6f Bohr") % _trust_radius).str() << flush;
            bool converged = false;

            // in 1st iteration, get the energy of the initial configuration
            if (_iteration == 0) {
                _last_energy = GetEnergy(); // in Hartree
                WriteTrajectory();
            }

            while (!converged) {
                _iteration++;
                _update_hessian = true; // always update Hessian after each accepted step
                _force_engine.Calculate(_last_energy);
                _force = _force_engine.GetForces();
                bool step_accepted = false;
                while (!step_accepted) {
                    BFGSStep();
                    Vector2QMAtoms();
                    _new_energy = GetEnergy();
                    _energy_delta = _new_energy - _last_energy;
                    step_accepted=AcceptRejectStep();
                } 

                // after the step is accepted, we can shift the stored data
                _last_energy = _new_energy;
                _old_xyz = _current_xyz;
                _current_xyz = _trial_xyz;
                _force_old = _force;

                WriteTrajectory();
                converged = CheckConvergence();
                Report();

                if (_iteration == _max_iteration) {
                    CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: not converged after %2$d iterations ") % _iteration % _max_iteration).str() << flush;
                    break;
                }

            } // convergence
            return;
        }

        /* Accept/reject the new geometry and adjust trust radius, if required */
        bool BFGSTRM::AcceptRejectStep() {
            bool step_accepted = false;
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
                step_accepted = true;
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: step accepted ") % _iteration).str() << flush;
                // adjust trust radius, if required
                double tr_check = _energy_delta / QuadraticEnergy();
                double norm_delta_pos=_delta_pos.norm();
                if (tr_check > 0.75 && 1.25 * norm_delta_pos > _trust_radius) {
                    _trust_radius = 2.0 * _trust_radius;
                    CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.6f Bohr") % _iteration % _trust_radius).str() << flush;
                } else if (tr_check < 0.25) {
                    _trust_radius = 0.25 * norm_delta_pos;
                    CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.6f Bohr") % _iteration % _trust_radius).str() << flush;
                }
            }
            return step_accepted;
        }

        /* Get the energy for the current configuration */
        double BFGSTRM::GetEnergy() {
            _gwbse_engine.setRedirectLogger(true);
            string logger_file = "gwbse_iteration_" + (boost::format("%1%") % _iteration).str() + ".log";
            _gwbse_engine.setLoggerFile(logger_file);
            _gwbse_engine.ExcitationEnergies(_orbitals);
            double energy = _orbitals.getTotalStateEnergy(_filter.CalcStateAndUpdate(_orbitals)); // in Hartree
            _gwbse_engine.setRedirectLogger(false);
            return energy;
        }
/*
        //Report results of accepted step
        void BFGSTRM::Report() {

            // accepted step
            CTP_LOG(ctp::logINFO, *_pLog) << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format(" =========== OPTIMIZATION SUMMARY ================================= ")).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << " At iteration  " << _iteration << flush;
            _force_engine.Report();

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   ---- POSITIONS (Angstrom)   ")).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Atom\t x\t  y\t  z ")).str() << flush;
            for (unsigned i = 0; i < ; i++) {
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format(" %1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f")
                        % i % (_current_xyz(i, 0) * votca::tools::conv::bohr2ang) % (_current_xyz(i, 1) * votca::tools::conv::bohr2ang) % (_current_xyz(i, 2) * votca::tools::conv::bohr2ang)).str() << flush;
            }

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Total energy:     %1$12.8f Hartree ") % _new_energy).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   quadratic energy: %1$12.8f Hartree ") % QuadraticEnergy()).str() << flush;
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
 */

        string BFGSTRM::Converged(bool converged) {
            if (converged) {
                return "converged    ";
            } else {
                return "not converged";
            }
        }

        bool BFGSTRM::CheckConvergence() {

            _energy_converged = false;
            _RMSForce_converged = false;
            _MaxForce_converged = false;
            _RMSStep_converged = false;
            _MaxStep_converged = false;

            Eigen::VectorXd xyz_shift = _current_pos - _previous_pos;

            _RMSForce = _current_gradient.cwiseAbs2().sum()/_current_gradient.size();
            _MaxForce = _current_gradient.cwiseAbs().maxCoeff();
            _RMSStep = xyz_shift.cwiseAbs2().sum()/xyz_shift.size();
            _MaxStep = xyz_shift.cwiseAbs().maxCoeff();

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



        /* Determine new trial coordinates according to BFGS */
        void BFGSTRM::BFGSStep() {
            WriteMatrixToVector();
            if (_update_hessian) UpdateHessian();
            PredictDisplacement();
            if (_delta_pos.norm()>_trust_radius){
                RegularizeStep();
            }
            QuadraticEnergy();
            WriteCoordinates2Matrices();
            return;
        }

      

        void BFGSTRM::UpdateHessian() {
            _delta_pos = _current_pos - _previous_pos;
            // we have no Hessian in the first iteration => start with something
            if (_iteration == 1) {
             _hessian=Eigen::MatrixXd::Identity(_delta_pos.size(),_delta_pos.size());
            }else {
                /* for later iteration, we can make use of an iterative refinement of
                 * the initial Hessian based on the gradient (force) history
                 */
                Eigen::VectorXd delta_gradient = _current_gradient - _previous_gradient;
                // second term in BFGS update (needs current Hessian)
               _hessian -= _hessian*_delta_pos*_delta_pos.transpose()*_hessian.transpose() / (_delta_pos.transpose()*_hessian*_delta_pos).value();
                // first term in BFGS update
                _hessian += (delta_gradient* delta_gradient.transpose()) / (delta_gradient.transpose()*_delta_pos);
                // symmetrize Hessian (since d2E/dxidxj should be symmetric)
               _hessian=0.5*(_hessian+_hessian.transpose());
            }
            return;
        }

        /* Predict displacement of atom coordinates */
        Eigen::VectorXd BFGSTRM::PredictDisplacement(const Eigen::MatrixXd& hessian,const Eigen::MatrixXd& gradient) const{
            return hessian.colPivHouseholderQr().solve(-gradient);
        }

        /* Check if predicted displacement leaves trust region */
        bool BFGSTRM::OutsideTrustRegion(double step) const{
            double trust_radius_squared = _trust_radius * _trust_radius;
            if (step > trust_radius_squared) {
                return true;
            } else {
                return false;
            }
        }

        /* Regularize step in case of prediction outside of Trust Region */
        Eigen::VectorXd BFGSTRM::RegularizeStep(const Eigen::VectorXd& delta_pos,const Eigen::VectorXd& gradient) const{

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_hessian);

            // start value for lambda  a bit lower than lowest eigenvalue of Hessian
            double lambda= 1.05 * es.eigenvalues()(0);
            if (es.eigenvalues()(0) > 0.0) {
                lambda = -0.05 * std::abs(es.eigenvalues()(0));
            } 
            // for constrained step, we expect
            double max_step_squared =delta_pos.squaredNorm();
            while (OutsideTrustRegion(max_step_squared)) {
                max_step_squared = 0.0;
                lambda -= 0.05 * std::abs(es.eigenvalues()(0));
                for (unsigned i = 0; i < delta_pos.size(); i++) {
                    double temp = es.eigenvectors().col(i).transpose()* gradient;
                    max_step_squared += temp * temp / (es.eigenvalues()(i) - lambda) / (es.eigenvalues()(i) - lambda);
                }
            }
                
            Eigen::VectorXd new_delta_pos = Eigen::VectorXd::Zero(delta_pos.size());
            for (unsigned i = 0; i < delta_pos.size(); _i++) {
                new_delta_pos -= es.eigenvectors().col(i) * (es.eigenvectors().col(i).transpose()*gradient) / (es.eigenvalues()(i) - lambda);
            }
            return new_delta_pos;
        }

        /* Estimate energy change based on quadratic approximation */
        double BFGSTRM::QuadraticEnergy() const{
            return (_current_gradient.transpose()* _delta_pos).value() + (0.5 * _delta_pos.transpose()*_hessian*_delta_pos).value();
        }
        
        
        Eigen::VectorXd BFGSTRM::Write3XMatrixToVector(const Eigen::MatrixX3d& matrix) const{
        Eigen::VectorXd vec=Eigen::VectorXd::Zero(matrix.rows()*matrix.cols());
        for (int i_cart = 0; i_cart < 3; i_cart++) {
            for (int i_atom = 0; i_atom < matrix.cols(); i_atom++) {
                    int idx = 3 * i_atom + i_cart;
                    vec(idx)=matrix(i_atom,i_cart);
                }
            }
            return vec;
        }

        Eigen::VectorXd BFGSTRM::QMAtomsToVector(std::vector<QMAtoms*>& atoms) const{
            Eigen::VectorXd pos=Eigen::VectorXd::Zero(3*atoms.size());
           for (unsigned i=0;i<atoms.size();i++) {
               pos(3*i) = atoms[i]->getPos().getX();
               pos(3*i+1)= atoms[i]->getPos().getY();
               pos(3*i+2)= atoms[i]->getPos().getZ();
            }
            return pos;
        }
        
        void BFGSTRM::Vector2QMAtoms(const Eigen::VectorXd& pos,std::vector<QMAtoms*>& atoms) const{
            for ( unsigned i = 0; i< atoms.size(); i++) {
                tools::vec pos_displaced(pos(3*i), pos(3*i+1) , pos(3*i+2));
                atoms[i]->setPos(pos_displaced);
            }
            return;
        }

        void BFGSTRM::WriteTrajectory(const std::string& filename,std::vector< QMAtom* >& atoms)const{
            ofstream ofs;
            if (_iteration == 0) {
                ofs.open(filename.c_str(), ofstream::out);
            } else {
                ofs.open(filename.c_str(), ofstream::app);
            }
            if (!ofs.is_open()) {
                throw runtime_error("Bad file handle: " + filename);
            }
            // write coordinates as xyz file
            ofs << atoms.size() << endl;
            ofs << "iteration " << _iteration << " energy " << _last_energy << " Hartree" << endl;

            for (const QMAtom* atom:atoms) {
                tools::vec pos=atom->getPos()*tools::conv::bohr2ang;
                ofs << atom->getType() << " " << pos.getX() << " " << pos.getY() << " " << pos.getZ() << endl;
            }
            ofs.close();
            return;
        }

    }
}
