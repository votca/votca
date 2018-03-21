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
#include "votca/xtp/convergenceacc.h"




namespace votca { namespace xtp {
   namespace ub = boost::numeric::ublas;
   
   
    Eigen::MatrixXd ConvergenceAcc::Iterate(const Eigen::MatrixXd& dmat,const Eigen::MatrixXd& H,Eigen::VectorXd &MOenergies,Eigen::MatrixXd &MOs,double totE){
      Eigen::MatrixXd H_guess=Eigen::MatrixXd::Zero(H.rows(),H.cols());    
    
      if(_mathist.size()>_histlength){
          delete _mathist[_maxerrorindex];
          delete _dmatHist[_maxerrorindex];
               _totE.erase(_totE.begin()+_maxerrorindex);
              _mathist.erase(_mathist.begin()+_maxerrorindex);
              _dmatHist.erase(_dmatHist.begin()+_maxerrorindex);
              
          }
          
      _totE.push_back(totE);
      
      Eigen::MatrixXd errormatrix=(*Sminusahalf).transpose()*(H*dmat*(*S)-(*S)*dmat*H)*(*Sminusahalf);
      _diiserror=errormatrix.cwiseAbs().maxCoeff();
      if(_noisy){
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " DIIs error " << _diiserror << flush;
            }
      Eigen::MatrixXd* old=new Eigen::MatrixXd;     
      *old=H;         
       _mathist.push_back(old);     
        Eigen::MatrixXd* dold=new Eigen::MatrixXd;     
      *dold=dmat;         
       _dmatHist.push_back(dold);
             
       if(_maxout){
        if (_diiserror>_maxerror){
            _maxerror=_diiserror;
            _maxerrorindex=_mathist.size();
        }
      } 
       
    diis.Update(_maxerrorindex,errormatrix);
    bool diis_error=false; 
   
       
    if (_diiserror<_adiis_start && _usediis && _mathist.size()>3){
        Eigen::VectorXd coeffs;
        //use ADIIs if energy has risen a lot in current iteration

        if(_diiserror>_diis_start || _totE[_totE.size()-1]>0.9*_totE[_totE.size()-2]){
            coeffs=adiis.CalcCoeff(_dmatHist,_mathist);
            diis_error=!adiis.Info();
            if(_noisy){
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using ADIIS" << flush;
            }
        }
        else if(_diiserror>0.0001 && _diiserror<_diis_start){
            Eigen::VectorXd coeffs1=diis.CalcCoeff();
            Eigen::VectorXd coeffs2=adiis.CalcCoeff(_dmatHist,_mathist);
            diis_error=!(diis.Info() && adiis.Info());
            double mixing=_diiserror/_diis_start;
            coeffs=mixing*coeffs2+(1-mixing)*coeffs1;
            if(_noisy){
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using ADIIS+DIIS" << flush;
            }
        }
        else{
             coeffs=diis.CalcCoeff();
             diis_error=!diis.Info();
             if(_noisy){
             CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using DIIS" << flush;
             }
        }
        if(diis_error){
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " DIIS failed using mixing instead" << flush;
          H_guess=H;
        }else{
        for (unsigned i=0;i<coeffs.size();i++){  
            if(std::abs(coeffs(i))<1e-8){ continue;}
            H_guess+=coeffs(i)*(*_mathist[i]);
            }
        }
           
    }
    else{       
        H_guess=H;     
    }
      
    double gap=MOenergies(_nocclevels)-MOenergies(_nocclevels-1);
      
    if((_diiserror>_levelshiftend && _levelshift>0.0) || gap<1e-6){
      if(_mode!=KSmode::fractional){
        Levelshift(H_guess,MOs);
      }
    }
    SolveFockmatrix( MOenergies,MOs,H_guess);
    Eigen::MatrixXd dmatout=DensityMatrix(MOs,MOenergies);
    
    mix.Updatemix(dmat,dmatout);
    if(_diiserror>_adiis_start ||!_usediis || diis_error ||_mathist.size()<=3 ){
      dmatout=mix.MixDmat(dmat,dmatout);
      if(_noisy){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using Mixing with alpha="<<mix.getAlpha() << flush;
        }
    }
    return dmatout;
    }
      
    void ConvergenceAcc::SolveFockmatrix(Eigen::VectorXd& MOenergies,Eigen::MatrixXd& MOs,const Eigen::MatrixXd&H){
        //transform to orthogonal for
        Eigen::MatrixXd H_ortho=(*Sminusahalf).transpose()*H*(*Sminusahalf);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_ortho);
        std::cout<<(es.info()==Eigen::ComputationInfo::Success)<<std::endl;
        
        MOenergies=es.eigenvalues();
        MOs=(*Sminusahalf)*es.eigenvectors();
        return;
    }
    
    void ConvergenceAcc::Levelshift(Eigen::MatrixXd& H,const Eigen::MatrixXd&MOs) {
        
        Eigen::MatrixXd virt = Eigen::MatrixXd::Zero(MOs.rows(),MOs.cols());
        for (unsigned _i = _nocclevels; _i < H.rows(); _i++) {
                        virt(_i, _i) = _levelshift; 
            }
        
        Eigen::PartialPivLU<Eigen::MatrixXd> lu(MOs);
        Eigen::MatrixXd inv=lu.inverse();
        if(_noisy){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using levelshift:" << _levelshift << " Hartree" << flush;
        }
        Eigen::MatrixXd vir=  inv.transpose()*virt*inv ; 
        H+=vir;
      
          return;
    }
    
    
    Eigen::MatrixXd ConvergenceAcc::DensityMatrix(const Eigen::MatrixXd& MOs, const Eigen::VectorXd& MOEnergies){
      Eigen::MatrixXd result;
      if(_mode== KSmode::closed){
        result=DensityMatrixGroundState( MOs);
      }else if(_mode==KSmode::open){
        result=DensityMatrixGroundState_unres( MOs);
      }else if(_mode==KSmode::fractional){
        result=DensityMatrixGroundState_frac(MOs, MOEnergies);
      }
      return result;
    }
    
    Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState(const Eigen::MatrixXd& MOs){
          const Eigen::MatrixXd occstates=MOs.block(0,0,MOs.rows(),_nocclevels);
          Eigen::MatrixXd dmatGS = 2.0*occstates*occstates.transpose();
          return dmatGS;
        } 
    
    
    Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState_unres(const Eigen::MatrixXd& MOs) {
            if (_nocclevels == 0) {
                return Eigen::MatrixXd::Zero(MOs.cols(),MOs.rows());
            }
            Eigen::MatrixXd occstates=MOs.block(0,0,MOs.rows(),_nocclevels);
            Eigen::MatrixXd dmatGS = occstates*occstates.transpose();
              return dmatGS;
        }

        Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState_frac(const Eigen::MatrixXd& MOs, const Eigen::VectorXd& MOEnergies) {
            if (_nocclevels == 0) {
                return Eigen::MatrixXd::Zero(MOs.rows(),MOs.cols());
            }
            unsigned numofelec=2*_nocclevels;
            Eigen::VectorXd occupation = Eigen::VectorXd::Zero(MOEnergies.size());

            double buffer = 0.0001;
            double homo_energy = MOEnergies(numofelec - 1);
            std::vector<unsigned> degeneracies;

            for (unsigned _level = 0; _level < occupation.size(); _level++) {
                if (MOEnergies(_level)<(homo_energy - buffer)) {
                    occupation(_level) = 1.0;
                    numofelec--;
                } else if (std::abs(MOEnergies(_level) - homo_energy) < buffer) {
                    degeneracies.push_back(_level);
                } else if (MOEnergies(_level)>(homo_energy + buffer)) {
                    occupation(_level) = 0.0;
                }
            }
            double deg_occupation = double(numofelec) / double(degeneracies.size());
            for (unsigned _level = 0; _level < degeneracies.size(); _level++) {
                occupation(degeneracies[_level]) = deg_occupation;
            }
            Eigen::MatrixXd _dmatGS =MOs*occupation.asDiagonal()*MOs.transpose();
            return _dmatGS;
        }
    
    
    
      
  
      
          
      

}}
