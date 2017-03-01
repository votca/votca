/*
 *            Copyright 2009-2016 The VOTCA Development Team
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
#include "votca/xtp/diis.h"

#if defined(GSL)
    #include "gsl/gsl_multimin.h"
#endif


namespace votca { namespace xtp {
   namespace ub = boost::numeric::ublas;
   
   
    double Diis::Evolve(const ub::matrix<double>& dmat,const ub::matrix<double>& H,ub::vector<double> &MOenergies,ub::matrix<double> &MOs, int this_iter, double totE){
      ub::matrix<double>H_guess=ub::zero_matrix<double>(H.size1(),H.size2());    
      if(_errormatrixhist.size()>_histlength){
          delete _mathist[_maxerrorindex];
          delete _errormatrixhist[_maxerrorindex];
          delete _Diis_Bs[_maxerrorindex];
          delete _dmathist[_maxerrorindex];
                _totE.erase(_totE.begin()+_maxerrorindex);
              _mathist.erase(_mathist.begin()+_maxerrorindex);
              _dmathist.erase(_mathist.begin()+_maxerrorindex);
              _errormatrixhist.erase(_errormatrixhist.begin()+_maxerrorindex);
              _Diis_Bs.erase( _Diis_Bs.begin()+_maxerrorindex);
              for( std::vector< std::vector<double>* >::iterator it=_Diis_Bs.begin();it<_Diis_Bs.end();++it){
                  std::vector<double>* vect=(*it);
                  vect->erase(vect->begin()+_maxerrorindex);
              }
               _FDs.erase(_FDs.begin()+_maxerrorindex);
              for( std::vector< std::vector<double>* >::iterator it=_FDs.begin();it<_FDs.end();++it){
                  std::vector<double>* vect=(*it);
                  vect->erase(vect->begin()+_maxerrorindex);
              }
          
          }
          
      _totE.push_back(totE);
   
     
      //cout<<"MOs"<< _orbitals->MOCoefficients()<< endl;
     
      //cout<<"D\n"<<_dftAOdmat<<endl; 
      //Calculate errormatrix and orthogonalize
      ub::matrix<double>temp=ub::prod(H,dmat);
      
      //cout<<"S\n"<<_dftAOoverlap._aomatrix<<endl; 
      ub::matrix<double> errormatrix=ub::prod(temp,(*S));
      //cout <<"FDS"<<endl;
      //cout << errormatrix<<endl;
      temp=ub::prod(dmat,H);
      //cout <<"SDF"<<endl;
      //cout<<ub::prod(_dftAOoverlap._aomatrix,temp)<<endl; 
      errormatrix-=ub::prod((*S),temp);
       //cout<<"before"<<endl;
     //cout<<errormatrix<<endl;
      temp=ub::prod(errormatrix,*Sminusahalf);
      errormatrix=ub::prod(ub::trans(*Sminusahalf),temp);
      
      temp.resize(0,0);
      //cout<<"after"<<endl;
      //cout<<errormatrix<<endl;
      
      double max=linalg_getMax(errormatrix,true);
   
      ub::matrix<double>* old=new ub::matrix<double>;     
      //exit(0);
      
      *old=H;         
       _mathist.push_back(old);
       
        ub::matrix<double>* dold=new ub::matrix<double>;     
      //exit(0);
      
      *dold=dmat;         
       _dmathist.push_back(dold);
  
      ub::matrix<double>* olderror=new ub::matrix<double>; 
      *olderror=errormatrix;
       _errormatrixhist.push_back(olderror);
       if(_maxout){
         
          if (max>_maxerror){
              _maxerror=max;
              _maxerrorindex=_mathist.size();
          }
      } 
       
      std::vector<double>* Bijs=new std::vector<double>;
       _Diis_Bs.push_back(Bijs);
      for (unsigned i=0;i<_errormatrixhist.size()-1;i++){
          double value=linalg_traceofProd(errormatrix,ub::trans(*_errormatrixhist[i]));
          Bijs->push_back(value);
          _Diis_Bs[i]->push_back(value);
      }
      Bijs->push_back(linalg_traceofProd(errormatrix,ub::trans(errormatrix)));
       
      if (max<_diis_start && _usediis && this_iter>2){
          bool _useold=false;
          bool check=false;
          //old Pulat DIIs
          
          ub::vector<double> coeffs=ub::zero_vector<double>(_mathist.size());
          if(_useold) {
          
          ub::matrix<double> B=ub::zero_matrix<double>(_mathist.size()+1);
          ub::vector<double> a=ub::zero_vector<double>(_mathist.size()+1);
          a(0)=-1;
          for (unsigned i=1;i<B.size1();i++){
              B(i,0)=-1;
              B(0,i)=-1;
          }
          //cout <<"Hello"<<endl;
          //cout<<"_errormatrixhist "<<_errormatrixhist.size()<<endl;
          //#pragma omp parallel for
          for (unsigned i=1;i<B.size1();i++){
              for (unsigned j=1;j<=i;j++){
                  //cout<<"i "<<i<<" j "<<j<<endl;
                  B(i,j)=_Diis_Bs[i-1]->at(j-1);
                  if(i!=j){
                    B(j,i)=B(i,j);
                  }
              }
          }
          //cout <<"solve"<<endl;
          check=linalg_solve(B,a);
          for(unsigned i=0;i<coeffs.size();i++){
              coeffs(i)=a(i+1);
          }
          }
          else{
          
              // C2-DIIS
          
            
            ub::matrix<double> B=ub::zero_matrix<double>(_mathist.size());
            ub::vector<double> eigenvalues=ub::zero_vector<double>(_mathist.size());
            ub::matrix<double> eigenvectors=ub::zero_matrix<double>(_mathist.size());
            
          
          for (unsigned i=0;i<B.size1();i++){
              for (unsigned j=0;j<=i;j++){
                  //cout<<"i "<<i<<" j "<<j<<endl;
                  B(i,j)=_Diis_Bs[i]->at(j);
                  if(i!=j){
                    B(j,i)=B(i,j);
                  }
              }
          }
          //cout<<"B:"<<B<<endl;
          linalg_eigenvalues(B,eigenvalues,eigenvectors);
          //cout<<"eigenvectors_b:"<<eigenvectors<<endl;
          // Normalize weights
          
        
          
         for (unsigned i=0;i<B.size1();i++){
         double norm=0.0;

         for (unsigned j=0;j<B.size2();j++){
         norm+=eigenvectors(j,i);    
        
         }
       
         for (unsigned j=0;j<B.size2();j++){
            eigenvectors(j,i)=eigenvectors(j,i)/norm;    
         }
          }
          //cout<<"eigenvectors_a:"<<eigenvectors<<endl;
          // Choose solution by picking out solution with smallest error
          ub::vector<double> errors=ub::zero_vector<double>(B.size1());
          ub::matrix<double> temp=ub::prod(B,eigenvectors);
          ub::matrix<double> eq=ub::prod(ub::trans(eigenvectors),temp);
          for (unsigned i=0;i<eq.size1();i++){
      
            errors(i)=eq(i,i);
          }
          //cout<<"Errors:"<<eq<<endl;
          double MaxWeight=10.0;
          //cout<<"eigenvectors"<<eigenvectors<<endl;
          double min=std::numeric_limits<double>::max();
          int minloc=-1;
          
            for (int i = 0; i < errors.size(); i++) {
                if (std::abs(errors(i)) < min) {

                    bool ok = true;
                    for (unsigned k = 0; k < eigenvectors.size2(); k++) {
                        if (eigenvectors(k, i) > MaxWeight) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        min = std::abs(errors(i));
                        minloc = int(i);
                    }
                }
            }
    
      if(minloc!=-1){
          check=true;
     for(unsigned k=0;k<eigenvectors.size2();k++){
       coeffs(k)=eigenvectors(k,minloc);   
     }
      }
      else{
          check=false;
       }
          }
     
          
          if (!check){
              H_guess=0.5*(*_mathist[0])+0.5*(*_mathist[1]);
              LOG(logDEBUG, *_pLog) << TimeStamp() << " Solving DIIs failed, just use mixing " << flush;
               
          }
          else{
                
                 for (unsigned i=0;i<_mathist.size();i++){  
                     if(std::abs(coeffs(i))<1e-8){ continue;}
                    H_guess+=coeffs(i)*(*_mathist[i]);
                    //cout <<i<<" "<<a(i+1,0)<<" "<<(*_mathist[i])<<endl;
                 }
                //cout <<"H_guess"<<H_guess<<endl;
          }
        
              
         
      }
      else{       
          H_guess=H;     
          }
      
      double gap=MOenergies(_nocclevels)-MOenergies(_nocclevels-1);
      
      if(max>0.2 || (_levelshift>0.00001 && gap<0.01)){
        Levelshift(H_guess,dmat,_levelshift,_unrestricted);
      }
      SolveFockmatrix( MOenergies,MOs,H_guess);
      //cout<<"after"<<MOenergies<<endl;
      //cout<<endl;
      return max;
      }
      
      void Diis::SolveFockmatrix(ub::vector<double>& MOenergies,ub::matrix<double>& MOs,const ub::matrix<double>&H){
          
           ub::matrix<double> temp;
           bool info=linalg_eigenvalues_general( H,(*S), MOenergies, temp);
            if (!info){
                    throw runtime_error("Generalized eigenvalue problem did not work.");
                }
           
           MOs=ub::trans(temp);
           
           return;
      }
      
      void Diis::Levelshift(ub::matrix<double>& H,const ub::matrix<double> & dmat,double levelshift,bool unrestricted){
          ub::matrix<double> I=ub::identity_matrix<double>(dmat.size1());
          double occupation=2.0;
          if(unrestricted){
              occupation=1.0;
          }
          H+=levelshift*(occupation*I-dmat);
                  
                  return;
      }

}}
