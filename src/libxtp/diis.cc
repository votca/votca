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
#include "votca/xtp/diis.h"




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
              _dmathist.erase(_dmathist.begin()+_maxerrorindex);
              _errormatrixhist.erase(_errormatrixhist.begin()+_maxerrorindex);
              _Diis_Bs.erase( _Diis_Bs.begin()+_maxerrorindex);
              for( std::vector< std::vector<double>* >::iterator it=_Diis_Bs.begin();it<_Diis_Bs.end();++it){
                  std::vector<double>* vect=(*it);
                  vect->erase(vect->begin()+_maxerrorindex);
              }
          }
          
      _totE.push_back(totE);

      ub::matrix<double>temp=ub::prod(H,dmat);
      ub::matrix<double> errormatrix=ub::prod(temp,(*S));
      temp=ub::prod(dmat,H);
      errormatrix-=ub::prod((*S),temp);
      temp=ub::prod(errormatrix,*Sminusahalf);
      errormatrix=ub::prod(ub::trans(*Sminusahalf),temp);
      
      temp.resize(0,0);
      
      double max=linalg_getMax(errormatrix,true);
   
      ub::matrix<double>* old=new ub::matrix<double>;     
      *old=H;         
       _mathist.push_back(old);
       
        ub::matrix<double>* dold=new ub::matrix<double>;     
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
         
      _DiF=ub::zero_vector<double>(_dmathist.size());
      _DiFj=ub::zero_matrix<double>(_dmathist.size());
       
    
  for(unsigned i=0;i<_dmathist.size();i++){
    _DiF(i)=linalg_traceofProd((*_dmathist[i])-dmat,H);
  }
  
  for(unsigned i=0;i<_dmathist.size();i++){
    for(unsigned j=0;j<_dmathist.size();j++){
        _DiFj(i,j)=linalg_traceofProd((*_dmathist[i])-dmat,(*_mathist[j])-H);
        }
   }
      
     
      
      
       
    if (max<_adiis_start && _usediis && this_iter>2){
        ub::vector<double> coeffs;
        //use EDIIs if energy has risen a lot in current iteration

        if(max>_diis_start || _totE[_totE.size()-1]>0.9*_totE[_totE.size()-2]){
            coeffs=ADIIsCoeff();
            if(_noisy){
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using ADIIS" << flush;
            }
        }
        else if(max>0.0001 && max<_diis_start){
            ub::vector<double> coeffs1=DIIsCoeff();
            //cout<<"DIIS "<<coeffs1<<endl;
            ub::vector<double> coeffs2=ADIIsCoeff();
            //cout<<"ADIIS "<<coeffs2<<endl;
            double mixing=max/_diis_start;
            coeffs=mixing*coeffs2+(1-mixing)*coeffs1;
            if(_noisy){
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using ADIIS+DIIS" << flush;
            }
        }
        else{
             coeffs=DIIsCoeff();
             if(_noisy){
             CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using DIIS" << flush;
             }
        }

       //check if last element completely rejected use mixing
        if(std::abs(coeffs(coeffs.size()-1))<0.001){ 
            coeffs=ub::zero_vector<double>(coeffs.size());
            coeffs(coeffs.size()-1)=0.3;
            coeffs(coeffs.size()-2)=0.7;
            if(_noisy){
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Last coefficient is too small use mixing with alpha=0.7 instead" << flush;
            }
            }


        for (unsigned i=0;i<coeffs.size();i++){  
            if(std::abs(coeffs(i))<1e-8){ continue;}
            H_guess+=coeffs(i)*(*_mathist[i]);
            //cout <<i<<" "<<a(i+1,0)<<" "<<(*_mathist[i])<<endl;
            }
            //cout <<"H_guess"<<H_guess<<endl;
    }
    else{       
        H_guess=H;     
    }
      
    double gap=MOenergies(_nocclevels)-MOenergies(_nocclevels-1);
      
    if((max>_levelshiftend && _levelshift>0.00001) || gap<1e-6){
        Levelshift(H_guess,MOs);
    }
    SolveFockmatrix( MOenergies,MOs,H_guess);
    return max;
    }
      
    void Diis::SolveFockmatrix(ub::vector<double>& MOenergies,ub::matrix<double>& MOs,ub::matrix<double>&H){
        //transform to orthogonal form
        
        ub::matrix<double>temp=ub::prod(ub::trans(*Sminusahalf),H);
        H=ub::prod(temp,*Sminusahalf);
        
        bool info=linalg_eigenvalues( MOenergies,H);
        if (!info){
            throw runtime_error("eigenvalue problem did not work.");
        }
       
        //H now stores the MOs
        MOs=ub::prod(*Sminusahalf,H);
        
        MOs=ub::trans(MOs);
        return;
    }
    
    void Diis::Levelshift(ub::matrix<double>& H,const ub::matrix<double>&MOs) {
        ub::matrix<double> transform=ub::trans(MOs);
        ub::matrix<double> trans_inv=ub::zero_matrix<double>(MOs.size1());
        linalg_invert(transform,trans_inv);
       
        ub::matrix<double> virt = ub::zero_matrix<double>(H.size1());
        for (unsigned _i = _nocclevels; _i < H.size1(); _i++) {
                        virt(_i, _i) = _levelshift; 
            }
        transform=ub::prod(ub::trans(trans_inv),virt);
        virt=ub::prod(transform,trans_inv);
        if(_noisy){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using levelshift:" << _levelshift << " Ha" << flush;
        }
        H +=  virt ; 
        
        return;
    }
      
      
    ub::vector<double>Diis::DIIsCoeff(){
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
          
            for (unsigned i = 0; i < errors.size(); i++) {
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
          
          if(!check){
               CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Solving DIIs failed, just use mixing " << flush;
               coeffs=ub::zero_vector<double>(_mathist.size());
               coeffs[coeffs.size()-1]=0.3;
               coeffs[coeffs.size()-2]=0.7;
          }
     return coeffs;  
   }   
      
      
      
   ub::vector<double> Diis::ADIIsCoeff(){
          
   size_t N=_DiF.size();
          
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf minfunc;
  minfunc.f = adiis::min_f;
  minfunc.df = adiis::min_df;
  minfunc.fdf = adiis::min_fdf;
  minfunc.n = N;
  minfunc.params = (void *) this;

  T=gsl_multimin_fdfminimizer_vector_bfgs2;
  s=gsl_multimin_fdfminimizer_alloc(T,N);

  // Starting point: equal weights on all matrices
    x=gsl_vector_alloc(N);
    gsl_vector_set_all(x,1.0/N);

  // Initialize the optimizer. Use initial step size 0.02, and an
  // orthogonality tolerance of 0.1 in the line searches (recommended
  // by GSL manual for bfgs).
  gsl_multimin_fdfminimizer_set(s, &minfunc, x, 0.02, 0.1);

  size_t iter=0;
  int status;
  do {
    iter++;
    //    printf("iteration %lu\n",iter);
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status) {
      //      printf("Error %i in minimization\n",status);
      break;
    }

    status = gsl_multimin_test_gradient(s->gradient, 1e-7);

    /*
    if (status == GSL_SUCCESS)
      printf ("Minimum found at:\n");
    printf("%5lu ", iter);
    for(size_t i=0;i<N;i++)
      printf("%.5g ",gsl_vector_get(s->x,i));
    printf("%10.5g\n",s->f);
    */
  }
  while (status == GSL_CONTINUE && iter < 1000);

  // Final estimate
  // double E_final=get_E(s->x);

  // Form minimum
  ub::vector<double> c=adiis::compute_c(s->x);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free (x);

  return c;
}

 
 
 
double Diis::get_E_adiis(const gsl_vector * x) const {
  // Consistency check
    if(x->size != _DiF.size()) {
        throw std::runtime_error("Incorrect number of parameters.");
    }

    ub::vector<double> c=adiis::compute_c(x);
    double Eval=0.0;
    for(unsigned i=0;i<c.size();i++){
        Eval+=2.0*c(i)*_DiF(i);
    }

    for(unsigned i=0;i<c.size();i++){
        for(unsigned j=0;j<c.size();j++){
            Eval+=c(i)*c(j)*_DiFj(i,j);
        }
    }

return Eval;
}

void Diis::get_dEdx_adiis(const gsl_vector * x, gsl_vector * dEdx) const {
  // Compute contraction coefficients
  ub::vector<double> c=adiis::compute_c(x);
  
   
  
  ub::vector<double> dEdc=2.0*_DiF + ub::prod(_DiFj,c) + ub::prod(ub::trans(_DiFj),c);
 

  // Compute jacobian of transformation: jac(i,j) = dc_i / dx_j
  ub::matrix<double> jac=adiis::compute_jac(x);

  // Finally, compute dEdx by plugging in Jacobian of transformation
  // dE/dx_i = dc_j/dx_i dE/dc_j
  ub::vector<double> dEdxv=ub::prod(ub::trans(jac),dEdc);
  for(size_t i=0;i< dEdxv.size();i++)
    gsl_vector_set(dEdx,i,dEdxv(i));
  return;
}

void Diis::get_E_dEdx_adiis(const gsl_vector * x, double * Eval, gsl_vector * dEdx) const {
  // Consistency check
   if(x->size != _DiF.size()) {
   
    throw std::runtime_error("Incorrect number of parameters.");
  }
  if(x->size != dEdx->size) {
    throw std::domain_error("x and dEdx have different sizes!\n");
  }

  // Compute energy
  *Eval=get_E_adiis(x);
  // and its derivative
  get_dEdx_adiis(x,dEdx);
  return;
}


ub::vector<double> adiis::compute_c(const gsl_vector * x) {
  // Compute contraction coefficients
  ub::vector<double> c=ub::zero_vector<double>(x->size);

  double xnorm=0.0;
  for(size_t i=0;i<x->size;i++) {
    c[i]=gsl_vector_get(x,i);
    c[i]=c[i]*c[i]; // c_i = x_i^2
    xnorm+=c[i];
  }
  for(size_t i=0;i<x->size;i++)
    c[i]/=xnorm; // c_i = x_i^2 / \sum_j x_j^2

  return c;
}  

ub::matrix<double> adiis::compute_jac(const gsl_vector * x) {
  // Compute jacobian of transformation: jac(i,j) = dc_i / dx_j

  // Compute coefficients
  std::vector<double> c(x->size);

  double xnorm=0.0;
  for(size_t i=0;i<x->size;i++) {
    c[i]=gsl_vector_get(x,i);
    c[i]=c[i]*c[i]; // c_i = x_i^2
    xnorm+=c[i];
  }
  for(size_t i=0;i<x->size;i++)
    c[i]/=xnorm; // c_i = x_i^2 / \sum_j x_j^2

 ub::matrix<double> jac=ub::zero_matrix<double>(c.size());
  for(size_t i=0;i<c.size();i++) {
    double xi=gsl_vector_get(x,i);

    for(size_t j=0;j<c.size();j++) {
      double xj=gsl_vector_get(x,j);

      jac(i,j)=-c[i]*2.0*xj/xnorm;
    }

    // Extra term on diagonal
    jac(i,i)+=2.0*xi/xnorm;
  }

  return jac;
}

double adiis::min_f(const gsl_vector * x, void * params) {
  Diis * a=(Diis *) params;
  return a->get_E_adiis(x);
}

void adiis::min_df(const gsl_vector * x, void * params, gsl_vector * g) {
  Diis * a=(Diis *) params;
  a->get_dEdx_adiis(x,g);
  return;
}

void adiis::min_fdf(const gsl_vector *x, void * params, double * f, gsl_vector * g) {
  Diis * a=(Diis *) params;
  a->get_E_dEdx_adiis(x,f,g);
  return;
}
      
      
          
      

}}
