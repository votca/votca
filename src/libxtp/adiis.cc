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
#include "votca/xtp/adiis.h"




namespace votca { namespace xtp {
   namespace ub = boost::numeric::ublas;
   
      
      
      
   Eigen::VectorXd ADIIS::CalcCoeff(const std::vector< Eigen::MatrixXd* >& _dmathist,const std::vector< Eigen::MatrixXd* >& _mathist){
      success=true;
      unsigned size=_dmathist.size();
      
      const Eigen::MatrixXd& dmat=*_dmathist[size-1];
      const Eigen::MatrixXd& H=*_mathist[size-1];
      _DiF = Eigen::VectorXd::Zero(size);
      _DiFj = Eigen::MatrixXd::Zero(size, size);


      for (unsigned i = 0; i < _dmathist.size(); i++) {
        _DiF(i) = ((*_dmathist[i]) - dmat).cwiseProduct(H).sum();
      }

      for (unsigned i = 0; i < _dmathist.size(); i++) {
        for (unsigned j = 0; j < _dmathist.size(); j++) {
          _DiFj(i, j) = ((*_dmathist[i]) - dmat).cwiseProduct((*_mathist[j]) - H).sum();
        }
      }   
     
     
std::cout<<_DiF<<std::endl;
std::cout<<_DiFj<<std::endl;     
          
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
  Eigen::VectorXd c=adiis::compute_c(s->x);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free (x);

  
  if(std::abs(c.tail(1).value())<0.001){     
        success=false;
      }
  return c;
}

 
 
 
double ADIIS::get_E_adiis(const gsl_vector * x) const {
  // Consistency check
    if(x->size != _DiF.size()) {
        throw std::runtime_error("Incorrect number of parameters.");
    }

    Eigen::VectorXd c=adiis::compute_c(x);
    double Eval=(2*c.transpose()*_DiF+c.transpose()*_DiFj*c).value();
    std::cout<<"E "<<Eval<<std::endl;

return Eval;
}

void ADIIS::get_dEdx_adiis(const gsl_vector * x, gsl_vector * dEdx) const {
  // Compute contraction coefficients
  Eigen::VectorXd c=adiis::compute_c(x);
  
   
  
  Eigen::VectorXd dEdc=2.0*_DiF + _DiFj*c + _DiFj.transpose()*c;
 

  // Compute jacobian of transformation: jac(i,j) = dc_i / dx_j
  Eigen::MatrixXd jac=adiis::compute_jac(x);

  // Finally, compute dEdx by plugging in Jacobian of transformation
  // dE/dx_i = dc_j/dx_i dE/dc_j
  Eigen::VectorXd dEdxv=jac.transpose()*dEdc;
  std::cout<<"Ed "<<dEdxv<<std::endl;
  for(size_t i=0;i< dEdxv.size();i++){
    gsl_vector_set(dEdx,i,dEdxv(i));
  }
  return;
}

void ADIIS::get_E_dEdx_adiis(const gsl_vector * x, double * Eval, gsl_vector * dEdx) const {
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


Eigen::VectorXd adiis::compute_c(const gsl_vector * x) {
  // Compute contraction coefficients
  Eigen::VectorXd c=Eigen::VectorXd::Zero(x->size);
  double xnorm=0.0;
  for(size_t i=0;i<x->size;i++) {
    c(i)=gsl_vector_get(x,i);
    xnorm+=c(i)*c(i);
    
  }
  c/=xnorm;
  return c;
}  

Eigen::MatrixXd adiis::compute_jac(const gsl_vector * x) {
  // Compute jacobian of transformation: jac(i,j) = dc_i / dx_j

  // Compute coefficients
  Eigen::VectorXd c=Eigen::VectorXd::Zero(x->size);

  double xnorm=0.0;
  for(size_t i=0;i<x->size;i++) {
    c(i)=gsl_vector_get(x,i);
    xnorm+=c(i)*c(i);
  }
  c/=xnorm;
  
 Eigen::MatrixXd jac=Eigen::MatrixXd::Zero(c.size(),c.size());
  for(size_t i=0;i<c.size();i++) {
    double xi=gsl_vector_get(x,i);

    for(size_t j=0;j<c.size();j++) {
      double xj=gsl_vector_get(x,j);

      jac(i,j)=-c(i)*2.0*xj/xnorm;
    }

    // Extra term on diagonal
    jac(i,i)+=2.0*xi/xnorm;
  }

  return jac;
}

double adiis::min_f(const gsl_vector * x, void * params) {
  ADIIS * a=(ADIIS *) params;
  return a->get_E_adiis(x);
}

void adiis::min_df(const gsl_vector * x, void * params, gsl_vector * g) {
  ADIIS * a=(ADIIS *) params;
  a->get_dEdx_adiis(x,g);
  return;
}

void adiis::min_fdf(const gsl_vector *x, void * params, double * f, gsl_vector * g) {
  ADIIS * a=(ADIIS *) params;
  a->get_E_dEdx_adiis(x,f,g);
  return;
}
      
      
          
      

}}
