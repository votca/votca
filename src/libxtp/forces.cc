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

#include <boost/format.hpp>

#include <votca/tools/elements.h>

#include <votca/xtp/atom.h>
#include <votca/xtp/forces.h>
#include "votca/xtp/statefilter.h"

namespace votca {
    namespace xtp {

      using std::flush;
        void Forces::Initialize(tools::Property &options) {

            std::vector<std::string> choices = {"forward", "central"};
            _force_method = options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(".method", choices);

            _noisy_output = options.ifExistsReturnElseReturnDefault<bool>(".noisy", false); 
           
            if ((_force_method == "forward") || (_force_method == "central")) {
                _displacement = options.ifExistsReturnElseReturnDefault<double>(".displacement", 0.001); // Angstrom
            }
            _displacement*=tools::conv::ang2bohr;

            // check for force removal options
            choices = {"total", "none"};
            std::string _force_removal = options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(".removal", choices);
            if (_force_removal == "total") _remove_total_force = true;

            _natoms = _orbitals.QMAtoms().size();
            _forces =Eigen::MatrixX3d::Zero(_natoms,3);

            return;
        }

        void Forces::Calculate(double energy) {

            TLogLevel ReportLevel = _pLog->getReportLevel(); // backup report level
            if ( ! _noisy_output ){
                _pLog->setReportLevel(logERROR); // go silent for force calculations
            }
            std::vector<QMAtom*>& atoms=_orbitals.QMAtoms();
            for (unsigned atom_index=0;atom_index<atoms.size();atom_index++) {
                if ( _noisy_output ){
                    XTP_LOG(logINFO, *_pLog) << "FORCES--DEBUG working on atom " << atom_index<< flush;
                }
                Eigen::Vector3d atom_force=Eigen::Vector3d::Zero();
                // Calculate Force on this atom
                if (_force_method == "forward") atom_force=NumForceForward(energy, atom_index);
                if (_force_method == "central") atom_force=NumForceCentral(energy, atom_index);
                _forces.row(atom_index)=atom_force.transpose();
            }
            _pLog->setReportLevel(ReportLevel); // 
            if (_remove_total_force) RemoveTotalForce();
            return;
        }

        void Forces::Report() const{

            XTP_LOG(logINFO, *_pLog) << (boost::format("   ---- FORCES (Hartree/Bohr)   ")).str() << flush;
            XTP_LOG(logINFO, *_pLog) << (boost::format("        %1$s differences   ") % _force_method).str() << flush;
            XTP_LOG(logINFO, *_pLog) << (boost::format("        displacement %1$1.4f Angstrom   ") % (_displacement*tools::conv::bohr2ang)).str() << flush;
            XTP_LOG(logINFO, *_pLog) << (boost::format("   Atom\t x\t  y\t  z ")).str() << flush;

            for (unsigned _i = 0; _i < _forces.rows(); _i++) {
                XTP_LOG(logINFO, *_pLog) << (boost::format(" %1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f")
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
                displacement_vec[i_cart]=_displacement;
                // update the coordinate
                tools::vec pos_displaced = current_pos + displacement_vec;
                atom->setPos(pos_displaced); 
                _gwbse_engine.ExcitationEnergies(_orbitals);
                double energy_displaced = _orbitals.getTotalStateEnergy(_filter.CalcState(_orbitals));
                force(i_cart) = (energy - energy_displaced) / _displacement;
                 atom->setPos(current_pos); // restore original coordinate into segment
            } // Cartesian directions
            return force;
        }

        /* Calculate forces on atoms numerically by central differences */
        Eigen::Vector3d Forces::NumForceCentral(double energy,int atom_index) {
            QMAtom* atom= _orbitals.QMAtoms()[atom_index];
            const tools::vec current_pos =atom->getPos();
            Eigen::Vector3d force=Eigen::Vector3d::Zero();
            for (unsigned i_cart = 0; i_cart < 3; i_cart++) {
                if ( _noisy_output ){
                    XTP_LOG(logINFO, *_pLog) << "FORCES--DEBUG           Cartesian component " << i_cart << flush;
                }
                tools::vec displacement_vec(0, 0, 0);          
                displacement_vec[i_cart]=_displacement;
                // update the coordinate
                tools::vec pos_displaced = current_pos + displacement_vec;
                atom->setPos(pos_displaced); 
                _gwbse_engine.ExcitationEnergies(_orbitals);
                double energy_displaced_plus = _orbitals.getTotalStateEnergy(_filter.CalcState(_orbitals));
                // update the coordinate
                pos_displaced = current_pos - displacement_vec;
                atom->setPos(pos_displaced); 
                _gwbse_engine.ExcitationEnergies(_orbitals);
                double energy_displaced_minus = _orbitals.getTotalStateEnergy(_filter.CalcState(_orbitals));
                force(i_cart) = 0.5 * (energy_displaced_minus - energy_displaced_plus) / _displacement;
                atom->setPos(current_pos); // restore original coordinate into orbital
            }
            return force;
        }

        void Forces::RemoveTotalForce() {
            Eigen::Vector3d total_force = _forces.colwise().sum();
            for (unsigned i_atom = 0; i_atom < _natoms; i_atom++) {
                _forces.row(i_atom)-=total_force/double(_natoms);
            }
            return;
        }

     
    }
}
