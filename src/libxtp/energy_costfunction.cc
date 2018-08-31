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


#include <votca/xtp/energy_costfunction.h>
#include <votca/xtp/qmatom.h>

namespace votca {
  namespace xtp {

    double Energy_costfunction::EvaluateCost(const Eigen::VectorXd& parameters){
      Vector2QMAtoms(parameters, _orbitals.QMAtoms());
      _gwbse_engine.setRedirectLogger(true);
      std::string logger_file = "gwbse_iteration_" + (boost::format("%1%") % _iteration).str() + ".log";
      _gwbse_engine.setLoggerFile(logger_file);
      _gwbse_engine.ExcitationEnergies(_orbitals);
      _gwbse_engine.setRedirectLogger(false);
      _energy= _orbitals.getTotalStateEnergy(_filter.CalcStateAndUpdate(_orbitals)); // in Hartree
      return _energy;
    }
            
    Eigen::VectorXd Energy_costfunction::EvaluateGradient(const Eigen::VectorXd& parameters){
      Vector2QMAtoms(parameters, _orbitals.QMAtoms());
      _force_engine.Calculate(_energy);
      Eigen::VectorXd gradient = Write3XMatrixToVector(-_force_engine.GetForces());
      return gradient;
    }
    
   

    bool Energy_costfunction::Converged(const Eigen::VectorXd& delta_parameters,
            double delta_cost, const Eigen::VectorXd& gradient) {
      _iteration++;
      bool energy_converged = false;
      bool RMSForce_converged = false;
      bool MaxForce_converged = false;
      bool RMSStep_converged = false;
      bool MaxStep_converged = false;
      _convval.deltaE = delta_cost;
      _convval.RMSForce = gradient.cwiseAbs2().sum() / gradient.size();
      _convval.MaxForce = gradient.cwiseAbs().maxCoeff();
      _convval.RMSStep = delta_parameters.cwiseAbs2().sum() / delta_parameters.size();
      _convval.MaxStep = delta_parameters.cwiseAbs().maxCoeff();

      if (std::abs(delta_cost) < _convpara.deltaE) energy_converged = true;
      if (_convval.RMSForce < _convpara.RMSForce) RMSForce_converged = true;
      if (_convval.MaxForce < _convpara.MaxForce) MaxForce_converged = true;
      if (_convval.RMSStep < _convpara.RMSStep) RMSStep_converged = true;
      if (_convval.MaxStep < _convpara.MaxStep) MaxStep_converged = true;

      if (energy_converged && RMSForce_converged && MaxForce_converged && RMSStep_converged && MaxStep_converged) {
        return true;
      }
      return false;
    }

    void Energy_costfunction::Vector2QMAtoms(const Eigen::VectorXd& pos, std::vector<QMAtom*>& atoms) {
      for (unsigned i = 0; i < atoms.size(); i++) {
        tools::vec pos_displaced(pos(3 * i), pos(3 * i + 1), pos(3 * i + 2));
        atoms[i]->setPos(pos_displaced);
      }
    }

    Eigen::VectorXd Energy_costfunction::QMAtoms2Vector(std::vector<QMAtom*>& atoms) {
      Eigen::VectorXd result = Eigen::VectorXd::Zero(3 * atoms.size());

      for (unsigned i = 0; i < atoms.size(); i++) {
        result(3 * i) = atoms[i]->getPos().getX();
        result(3 * i + 1) = atoms[i]->getPos().getY();
        result(3 * i + 2) = atoms[i]->getPos().getZ();
      }
      return result;
    }

    Eigen::VectorXd Energy_costfunction::Write3XMatrixToVector(const Eigen::MatrixX3d& matrix) {
      Eigen::VectorXd vec = Eigen::VectorXd::Zero(matrix.rows() * matrix.cols());
      for (int i_cart = 0; i_cart < 3; i_cart++) {
        for (int i_atom = 0; i_atom < matrix.cols(); i_atom++) {
          int idx = 3 * i_atom + i_cart;
          vec(idx) = matrix(i_atom, i_cart);
        }
      }
      return vec;
    }

  }
}
