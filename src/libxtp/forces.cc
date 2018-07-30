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

#include <votca/xtp/forces.h>
#include <boost/format.hpp>

namespace votca {
    namespace xtp {

      using std::flush;
        void Forces::Initialize(tools::Property &options) {

            // pre-check forces method
            std::vector<std::string> choices = {"forward", "central"};
            _force_method = options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(".method", choices);

            // output level
            _noisy_output = options.ifExistsReturnElseReturnDefault<bool>(".noisy", false); 
            
            
            // precaution in case we implement approx. analytic forces in the future
            if ((_force_method == "forward") || (_force_method == "central")) {
                _displacement = options.ifExistsReturnElseReturnDefault<double>(".displacement", 0.001); // Angstrom
            }

            // check for force removal options
            choices = {"total", "CoM", "none"};
            std::string _force_removal = options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(".removal", choices);
            if (_force_removal == "total") _remove_total_force = true;
            if (_force_removal == "CoM") _remove_CoM_force = true;

            _natoms = _orbitals.QMAtoms().size();
            _forces =Eigen::MatrixX3d::Zero(_natoms,3);


            return;

        }

        void Forces::Calculate(double energy) {

            ctp::TLogLevel ReportLevel = _pLog->getReportLevel(); // backup report level
            if ( ! _noisy_output ){
                _pLog->setReportLevel(ctp::logERROR); // go silent for force calculations
            }

            std::vector<QMAtom*>& atoms=_orbitals.QMAtoms();
            for (unsigned atom_index=0;atom_index<atoms.size();atom_index++) {
                QMAtom* atom=atoms[atom_index];
                if ( _noisy_output ){
                    CTP_LOG(ctp::logINFO, *_pLog) << "FORCES--DEBUG working on atom " << atom_index<< flush;
                }
                Eigen::Vector3d atom_force;
                // Calculate Force on this atom
                
                if (_force_method == "forward") atom_force=NumForceForward(energy, atom_index);
                if (_force_method == "central") atom_force=NumForceCentral(energy, atom_index);
                
                _forces.row(atom_index)=atom_force.transpose();
            }
            _pLog->setReportLevel(ReportLevel); // 
            if (_remove_total_force) RemoveTotalForce();
            //Report();

            return;
        }

        void Forces::Report() {

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   ---- FORCES (Hartree/Bohr)   ")).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("        %1$s differences   ") % _force_method).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("        displacement %1$1.4f Angstrom   ") % _displacement).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Atom\t x\t  y\t  z ")).str() << flush;

            for (unsigned _i = 0; _i < _forces.rows(); _i++) {
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format(" %1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f")
                        % _i % _forces(_i, 0) % _forces(_i, 1) % _forces(_i, 2)).str() << flush;
            }

            return;

        }

        /* Calculate forces on an atom numerically by forward differences */
        Eigen::Vector3d Forces::NumForceForward(double energy,int atom_index) {
            Eigen::Vector3d force=Eigen::Vector3d::Zero();
            // get this atoms's current coordinates
            QMAtom* atom= _orbitals.QMAtoms()[atom_index];
            const tools::vec current_pos =atom->getPos();

            for (unsigned i_cart = 0; i_cart < 3; i_cart++) {

                tools::vec displacement_vec(0, 0, 0);          
                displacement_vec[i_cart]=_displacement * tools::conv::ang2bohr; //  _displacement in Angstrom
               
                // update the coordinate
                tools::vec pos_displaced = current_pos + displacement_vec;

                atom->setPos(pos_displaced); 

                // run DFT and GW-BSE for this geometry
                _gwbse_engine.ExcitationEnergies(_qmpackage, _orbitals);

                // get total energy for this excited state
                double energy_displaced = _orbitals.getTotalExcitedStateEnergy(_spin_type, _opt_state);

                // calculate force and put into matrix
                force(i_cart) = (energy - energy_displaced) / (_displacement * votca::tools::conv::ang2bohr); // force a.u./a.u.
                 atom->setPos(current_pos); // restore original coordinate into segment
            } // Cartesian directions
            return force;
        }

        /* Calculate forces on atoms numerically by central differences */
        Eigen::Vector3d Forces::NumForceCentral(double energy,int atom_index) {


            QMAtom* atom= _orbitals.QMAtoms()[atom_index];
            const tools::vec current_pos =atom->getPos();
            Eigen::Vector3d force=Eigen::Vector3d::Zero();
            // go through all cartesian components
            for (unsigned i_cart = 0; i_cart < 3; i_cart++) {

                if ( _noisy_output ){
                    CTP_LOG(ctp::logINFO, *_pLog) << "FORCES--DEBUG           Cartesian component " << i_cart << flush;
                }
                
                tools::vec displacement_vec(0, 0, 0);          
                displacement_vec[i_cart]=_displacement * tools::conv::ang2bohr; //  _displacement in Angstrom

                // update the coordinate
                tools::vec pos_displaced = current_pos + displacement_vec;

                atom->setPos(pos_displaced); 
                // run DFT and GW-BSE for this geometry
                _gwbse_engine.ExcitationEnergies(_qmpackage, _orbitals);

                // get total energy for this excited state
                double energy_displaced_plus = _orbitals.getTotalExcitedStateEnergy(_spin_type, _opt_state);

                // get displacement vector in negative direction

                // update the coordinate
                pos_displaced = current_pos - displacement_vec;
                atom->setPos(pos_displaced); 
                // run DFT and GW-BSE for this geometry
                _gwbse_engine.ExcitationEnergies(_qmpackage,_orbitals);

                // get total energy for this excited state
                double energy_displaced_minus = _orbitals.getTotalExcitedStateEnergy(_spin_type, _opt_state);
                // calculate force and put into matrix
                force(i_cart) = 0.5 * (energy_displaced_minus - energy_displaced_plus) / (_displacement * votca::tools::conv::ang2bohr); // force a.u./a.u.
                atom->setPos(current_pos); // restore original coordinate into segment
            }

            return force;
        }

        /* Adjust forces so that sum of forces is zero */
        void Forces::RemoveTotalForce() {
            // total force on all atoms
            Eigen::Vector3d _total_force = TotalForce();
            // zero total force
            for (unsigned _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                _forces.row(_i_atom)-=_total_force/double(_natoms);
            }
            return;
        }

        /* Determine Total Force on all atoms */
         Eigen::Vector3d Forces::TotalForce() {
            return _forces.colwise().sum();
        }
    }
}
