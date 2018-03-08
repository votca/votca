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
#include "votca/xtp/mixing.h"



namespace votca { namespace xtp {
   namespace ub = boost::numeric::ublas;
   
     Eigen::MatrixXd Mixing::MixDmat(const Eigen::MatrixXd& dmatin,const Eigen::MatrixXd& dmatout,bool noisy ){
          double alpha=0.0;
          Updatemix(dmatin,dmatout);
          if(!_automaticmixing){   
              alpha=_mixingparameter;
          }
          else{
              
          if(_Pout.size()<2){
              alpha=0.7;             
          }
          else{
              
              Eigen::VectorXd nominator=*(_Pout[1])-*(_Pout[0]);
              Eigen::VectorXd denominator=nominator-(*(_Pin[1])-*(_Pin[0]));
              //cout<<endl;
              //cout<<_Pout.size()<<endl;
              //cout<<nominator<<endl;
              //cout<<denominator<<endl;
              double count=0.0;
              for(unsigned i=0;i<nominator.size();i++){
                  double a=nominator(i)/denominator(i);
                  if(a<0.2){
                      a=0.2;
                  }
                  else if(a>0.8){
                      a=0.8;
                  }
                  count+=a;
              }
              alpha=count/double(nominator.size());
             
              
          }
          }
          if(noisy){
          CTP_LOG(votca::ctp::logDEBUG, *_pLog) << votca::ctp::TimeStamp() << " Using Mixing with mixingparamter="<<alpha<< flush;
          }
          
          Eigen::MatrixXd dmatnew=alpha*dmatin+(1.0-alpha)*dmatout;
          return dmatnew;
      }
      
      void Mixing::Updatemix(const Eigen::MatrixXd& dmatin,const Eigen::MatrixXd& dmatout ){
          if(_Pin.size()>1){
              delete _Pin[0];
              _Pin.erase( _Pin.begin());
          }
          if(_Pout.size()>1){
              delete _Pout[0];
              _Pout.erase( _Pout.begin());
          }
          Eigen::VectorXd* mcharges=new Eigen::VectorXd;
          (*mcharges)=Mullikencharges(dmatin);
          _Pin.push_back(mcharges);
          mcharges=new Eigen::VectorXd;
          (*mcharges)=Mullikencharges(dmatout);
          _Pout.push_back(mcharges);
      }
      
      Eigen::VectorXd Mixing::Mullikencharges(const Eigen::MatrixXd& dmat){
          return (dmat*(*S)).diagonal();
      }

}}
