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

#include <votca/xtp/bfgs-trm.h>
#include <votca/xtp/atom.h>
#include <boost/format.hpp>

namespace votca {
    namespace xtp {

        void BFGSTRM::Optimize(const Eigen::VectorXd& initialparameters) {
           _parameters=initialparameters;

            Eigen::VectorXd gradient=Eigen::VectorXd::Zero(_parameters.size());
            Eigen::VectorXd delta_p_trial=Eigen::VectorXd::Zero(_parameters.size());
            
            double lastcost=_costfunction.EvaluateCost(_parameters);
            Eigen::VectorXd last_gradient=Eigen::VectorXd::Zero(_parameters.size());
            double delta_cost=0;
            
            for(_iteration=0;_iteration<_max_iteration;_iteration++){
               gradient=_costfunction.EvaluateGradient(_parameters);
                bool step_accepted = false;
                for (int i=0;i<100;i++){
                    delta_p_trial=CalculateInitialStep(gradient);
                    if (delta_p_trial.norm()>_trust_radius){
                       delta_p_trial=CalculateRegularizedStep(delta_p_trial,gradient);
                    }
                    double trialcost=_costfunction.EvaluateCost(_parameters+delta_p_trial);
                    delta_cost=trialcost-lastcost;
                    step_accepted=AcceptRejectStep(delta_p_trial,gradient,delta_cost);
                    if(step_accepted){
                      _cost=trialcost;
                      _parameters+=delta_p_trial;
                      break;
                    }
                    
                } 
                if(_iteration>0){
                  UpdateHessian(delta_p_trial,gradient-last_gradient);
                }
                lastcost=_cost;
                last_gradient=gradient;
                for(auto& func:_callbacks){
                 func();
                }
                if(_costfunction.Converged(delta_p_trial,delta_cost,gradient)){
                  break;
                } else if(_iteration == _max_iteration-1) {
                  _success=false;
                  if(_logging){
                    XTP_LOG(logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: not converged after %2$d iterations ")
                            % _iteration % _max_iteration).str() << std::flush;
                  }
                }

            }
            return;
        }

        /* Accept/reject the new geometry and adjust trust radius, if required */
        bool BFGSTRM::AcceptRejectStep(const Eigen::VectorXd& delta_p,const Eigen::VectorXd& gradient,double cost_delta) {
            bool step_accepted = false;
            if (cost_delta > 0.0) {
                // total energy has unexpectedly increased, half the trust radius
                _trust_radius = 0.25 * _trust_radius;
                if(_logging){
                  XTP_LOG(logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: step rejected ")
                          % _iteration).str() << std::flush;
                  XTP_LOG(logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.6f") 
                          % _iteration % _trust_radius).str() << std::flush;
                }
            } else {
                // total energy has decreased, we accept the step but might update the trust radius
                step_accepted = true;
                // adjust trust radius, if required
                double tr_check = cost_delta / QuadraticEnergy(gradient,delta_p);
                double norm_delta_p=delta_p.squaredNorm();
                if (tr_check > 0.75 && 1.25 * norm_delta_p > _trust_radius*_trust_radius) {
                    _trust_radius = 2.0 * _trust_radius;
                } else if (tr_check < 0.25) {
                    _trust_radius = 0.25 * _trust_radius;
                }
                if(_logging){
                  XTP_LOG(logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: step accepted ")
                          % _iteration).str() << std::flush;
                  XTP_LOG(logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.6f")
                          % _iteration % _trust_radius).str() << std::flush;
                }
            }
            return step_accepted;
        }

        void BFGSTRM::UpdateHessian(const Eigen::VectorXd& delta_pos,const Eigen::VectorXd& delta_gradient) {
                // second term in BFGS update (needs current Hessian)
               _hessian -= _hessian*delta_pos*delta_pos.transpose()*_hessian.transpose() / (delta_pos.transpose()*_hessian*delta_pos).value();
                // first term in BFGS update
                _hessian += (delta_gradient* delta_gradient.transpose()) / (delta_gradient.transpose()*delta_pos);
                // symmetrize Hessian (since d2E/dxidxj should be symmetric)
               _hessian=0.5*(_hessian+_hessian.transpose());
            return;
        }

        /* Predict displacement of atom coordinates */
        Eigen::VectorXd BFGSTRM::CalculateInitialStep(const Eigen::MatrixXd& gradient)const{
            return _hessian.colPivHouseholderQr().solve(-gradient);
        }

        /* Regularize step in case of prediction outside of Trust Region */
        Eigen::VectorXd BFGSTRM::CalculateRegularizedStep(const Eigen::VectorXd& delta_pos,const Eigen::VectorXd& gradient) const{
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_hessian);
            // start value for lambda  a bit lower than lowest eigenvalue of Hessian
            double lambda= 1.05 * es.eigenvalues()(0);
            if (es.eigenvalues()(0) > 0.0) {
              lambda = -0.05 * std::abs(es.eigenvalues()(0));
            } 
            // for constrained step, we expect
            double max_step_squared =delta_pos.squaredNorm();
            while (max_step_squared>(_trust_radius * _trust_radius)) {
              lambda -= 0.05 * std::abs(es.eigenvalues()(0));
              Eigen::VectorXd quotient=(es.eigenvalues().array()-lambda).cwiseAbs2();
              auto factor=(es.eigenvectors().transpose()*gradient).cwiseAbs2();
              max_step_squared = (factor.cwiseQuotient(quotient)).sum();
            }
                
            Eigen::VectorXd new_delta_pos = Eigen::VectorXd::Zero(delta_pos.size());
            for (unsigned i = 0; i < delta_pos.size(); i++) {
                new_delta_pos -= es.eigenvectors().col(i) * (es.eigenvectors().col(i).transpose()*gradient) / (es.eigenvalues()(i) - lambda);
            }
            return new_delta_pos;
        }

        /* Estimate energy change based on quadratic approximation */
        double BFGSTRM::QuadraticEnergy(const Eigen::VectorXd& gradient, const Eigen::VectorXd& delta_pos) const{
            return (gradient.transpose()* delta_pos).value() + 0.5*(delta_pos.transpose()*_hessian*delta_pos).value();
        }
        
        


    }
}
