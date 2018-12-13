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
      std::string logger_file;
      if(_iteration==0){
          logger_file="gwbse_iteration_" + (boost::format("%1%") % _iteration).str() + ".log";
          _iteration++;
      }else{
         logger_file="gwbse_iteration_" + (boost::format("%1%") % _iteration).str() + ".log";
      }
      _gwbse_engine.setLoggerFile(logger_file);
      _gwbse_engine.ExcitationEnergies(_orbitals);
      _gwbse_engine.setRedirectLogger(false);
      _energy= _orbitals.getTotalStateEnergy(_filter.CalcStateAndUpdate(_orbitals)); // in Hartree
      return _energy;
    }
            
    Eigen::VectorXd Energy_costfunction::EvaluateGradient(const Eigen::VectorXd& parameters){
      Vector2QMAtoms(parameters, _orbitals.QMAtoms());
      _force_engine.Calculate(_orbitals);
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
      Energy_costfunction::conv_paras convval;
      convval.deltaE = delta_cost;
      convval.RMSForce =  std::sqrt(gradient.cwiseAbs2().sum()) / gradient.size();
      convval.MaxForce = gradient.cwiseAbs().maxCoeff(&convval.maxforceindex);
      convval.RMSStep = std::sqrt(delta_parameters.cwiseAbs2().sum()) / delta_parameters.size();
      convval.MaxStep = delta_parameters.cwiseAbs().maxCoeff(&convval.maxstepindex);

      if (std::abs(convval.deltaE) < _convpara.deltaE) energy_converged = true;
      if (convval.RMSForce < _convpara.RMSForce) RMSForce_converged = true;
      if (convval.MaxForce < _convpara.MaxForce) MaxForce_converged = true;
      if (convval.RMSStep < _convpara.RMSStep) RMSStep_converged = true;
      if (convval.MaxStep < _convpara.MaxStep) MaxStep_converged = true;
      Report(convval);
      if (energy_converged && RMSForce_converged && MaxForce_converged && RMSStep_converged && MaxStep_converged) {
        return true;
      }
      return false;
    }
    
   void Energy_costfunction::Report(const Energy_costfunction::conv_paras& val){
      const Energy_costfunction::conv_paras& paras = getConvParas();
     XTP_LOG(logINFO, *_pLog) << (boost::format("   energy change:    %1$12.8f Hartree      %2$s") % val.deltaE % Converged(val.deltaE, paras.deltaE)).str() << std::flush;
      XTP_LOG(logINFO, *_pLog) << (boost::format("   RMS force:        %1$12.8f Hartree/Bohr %2$s") % val.RMSForce % Converged(val.RMSForce, paras.RMSForce)).str() << std::flush;
      XTP_LOG(logINFO, *_pLog) << (boost::format("   Max force:        %1$12.8f Hartree/Bohr %2$s") % val.MaxForce % Converged(val.MaxForce, paras.MaxForce)).str() << std::flush;
      XTP_LOG(logINFO, *_pLog) << (boost::format("   RMS step:         %1$12.8f Bohr         %2$s") % val.RMSStep % Converged(val.RMSStep, paras.RMSStep)).str() << std::flush;
      XTP_LOG(logINFO, *_pLog) << (boost::format("   Max step:         %1$12.8f Bohr         %2$s") % val.MaxStep % Converged(val.MaxStep, paras.MaxStep)).str() << std::flush;
      XTP_LOG(logINFO, *_pLog) << (boost::format(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")).str() << std::flush;
      XTP_LOG(logINFO, *_pLog) << std::flush;
   }
   
   std::string Energy_costfunction::Converged(double val, double limit){
      if (std::abs(val) < limit) {
        return (boost::format("Converged (%1%)") % limit).str();
      } else {
        return (boost::format("Not converged (%1%)") % limit).str();
      }
    }
    
    

    void Energy_costfunction::Vector2QMAtoms(const Eigen::VectorXd& pos,QMMolecule& atoms) {
      for (int i = 0; i < atoms.size(); i++) {
        Eigen::Vector3d pos_displaced;
        pos_displaced<<pos(3 * i), pos(3 * i + 1), pos(3 * i + 2);
        atoms[i].setPos(pos_displaced);
      }
    }

    Eigen::VectorXd Energy_costfunction::QMAtoms2Vector(QMMolecule& atoms) {
      Eigen::VectorXd result = Eigen::VectorXd::Zero(3 * atoms.size());

      for (int i = 0; i < atoms.size(); i++) {
        result(3 * i) = atoms[i].getPos().x();
        result(3 * i + 1) = atoms[i].getPos().y();
        result(3 * i + 2) = atoms[i].getPos().z();
      }
      return result;
    }

    Eigen::VectorXd Energy_costfunction::Write3XMatrixToVector(const Eigen::MatrixX3d& matrix) {
      Eigen::VectorXd vec = Eigen::VectorXd::Zero(matrix.rows() * matrix.cols());
      for (int i_cart = 0; i_cart < 3; i_cart++) {
        for (int i_atom = 0; i_atom < matrix.rows(); i_atom++) {
          int idx = 3 * i_atom + i_cart;
          vec(idx) = matrix(i_atom, i_cart);
        }
      }
      return vec;
    }

  }
}
