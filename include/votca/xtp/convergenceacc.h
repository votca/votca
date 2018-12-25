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

#ifndef VOTCA_XTP_CONVERGENCEACC_H
#define VOTCA_XTP_CONVERGENCEACC_H

#include <votca/xtp/basisset.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/adiis.h>
#include <votca/xtp/diis.h>
#include <votca/xtp/aomatrix.h>

namespace votca { namespace xtp {

 
class ConvergenceAcc{
    public:
    
    enum KSmode { closed, open, fractional };

    ConvergenceAcc() {_mode=KSmode::closed;
                       _usemixing=true;
                      _diiserror=std::numeric_limits<double>::max();
                      _maxerrorindex=0;
                      _maxerror=0.0;};
                      
   ~ConvergenceAcc() {
     for (std::vector< Eigen::MatrixXd* >::iterator it = _mathist.begin() ; it !=_mathist.end(); ++it){
         delete *it;
     }
     _mathist.clear();
      for (std::vector< Eigen::MatrixXd* >::iterator it = _dmatHist.begin() ; it !=_dmatHist.end(); ++it){
         delete *it;
     }
     _dmatHist.clear();
   }
   
   void Configure(KSmode mode,bool usediis,bool noisy, 
                    int histlength, bool maxout, double adiis_start,
                    double diis_start,double levelshift,double levelshiftend,
                    int numberofelectrons, double mixingparameter){
       _mode=mode;
       _usediis=usediis;
       _noisy=noisy;
       _histlength=histlength;
       _diis.setHistLength(histlength);
       _maxout=maxout;
       _adiis_start=adiis_start;
       _diis_start=diis_start;
       _levelshift=levelshift;
       _levelshiftend=levelshiftend;
       _mixingparameter=mixingparameter;
       _numberofelectrons=numberofelectrons;
       if(mode==KSmode::closed){
           _nocclevels=_numberofelectrons/2;
       }
       else if(mode==KSmode::open){
           _nocclevels=_numberofelectrons;
       }
       else if(mode==KSmode::fractional){
           _nocclevels=0;
       }
   }
   
   void setOverlap(AOOverlap* S, double etol);
   
   double getDIIsError(){return _diiserror;}
   
    bool getUseMixing(){return _usemixing;}
   
    void setLogger(Logger *pLog){_pLog=pLog;}
    Eigen::MatrixXd Iterate(const Eigen::MatrixXd& dmat,Eigen::MatrixXd& H,Eigen::VectorXd &MOenergies,Eigen::MatrixXd &MOs,double totE);
    void SolveFockmatrix(Eigen::VectorXd& MOenergies,Eigen::MatrixXd& MOs,const Eigen::MatrixXd&H);
    void Levelshift(Eigen::MatrixXd& H);

    Eigen::MatrixXd DensityMatrix(const Eigen::MatrixXd& MOs, const Eigen::VectorXd& MOEnergies);
   
 private:
     
    KSmode _mode;
    Eigen::MatrixXd DensityMatrixGroundState(const Eigen::MatrixXd& MOs);
    Eigen::MatrixXd DensityMatrixGroundState_unres(const Eigen::MatrixXd& MOs);
    Eigen::MatrixXd DensityMatrixGroundState_frac(const Eigen::MatrixXd& MOs, const Eigen::VectorXd& MOEnergies);
     
    bool                                _usemixing;
    Logger *                       _pLog;
    const AOOverlap* _S;
    
    bool                              _usediis;
    bool                              _noisy;
    int                          _histlength;
    bool                              _maxout;
    Eigen::MatrixXd                 Sminusahalf;
    Eigen::MatrixXd                 Sonehalf;
    Eigen::MatrixXd                 MOsinv;
    double                              _maxerror;
    double                              _diiserror;
    double                              _adiis_start;  
    double                              _diis_start;
    double                              _levelshiftend;
    int                            _maxerrorindex;
    std::vector< Eigen::MatrixXd* >   _mathist;
    std::vector< Eigen::MatrixXd* >   _dmatHist;
    double                              _mixingparameter;
    std::vector<double>                 _totE;
   
    int _numberofelectrons;
    int _nocclevels;
    double _levelshift;
    
    ADIIS _adiis;
    DIIS _diis;
   
    
    
    
    
  
 };
 
 
}}

#endif // VOTCA_XTP_CONVERGENCEACC_H	

