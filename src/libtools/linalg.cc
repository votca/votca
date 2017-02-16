/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include <votca/tools/linalg.h>

#include <iostream>
#include <sstream>
namespace votca { namespace tools {

using namespace std;
namespace ub = boost::numeric::ublas;
int linalg_invert_svd(ub::matrix<double> &A, ub::matrix<double> &V,double limitCN){
    int dimensions=0;
    
    if(A.size1()!=A.size2()){
        throw runtime_error("You cannot invert a non quadratic matrix!. linalg_invert_svd");
    }
     if(A.size1()==0){
        throw runtime_error("Matrix has size 0. linalg_invert_svd");
    }
    //make local copy because svd destroys A
    ub::matrix<double> work=A;
    ub::matrix<double>VT=ub::zero_matrix<double>(A.size1(),A.size2());
    ub::vector<double>S=ub::zero_vector<double>(A.size1());
    
    bool check=linalg_singular_value_decomposition(work, VT,S );
    if (check){
        throw runtime_error("Svd Calculation failed. linalg_invert_svd");
    }
    /*
    cout << "Hallo"<<endl;
    for(int i=0;i<work.size1();i++){
        for(int j=0;j<work.size2();j++){
            cout << "U["<<i<<":"<<j<<"]="<<work(i,j)<<endl;
    }
    }
   */
   

    //invert S
    ub::matrix<double>Sinverse=ub::zero_matrix<double>(S.size());
    for (unsigned i=0;i<S.size();i++){
        if (S(i)==0){
            break;
        }
        double CN=S(0)/S(i);
        //cout<< CN<<endl;
        //cout << limitCN << endl;
        if(CN>limitCN){
            break;
        }
        else{
            Sinverse(i,i)=1.0/S(i);
            dimensions++;
        }
    }
    
    
    ub::matrix<double> _temp=ub::prod(ub::trans(VT),Sinverse);
    V=ub::prod(_temp,ub::trans(work));
    
    return S.size()-dimensions;
}




double linalg_loewdin(ub::matrix<double> &J, ub::matrix<double> &S){
    if (J.size1()!=J.size2() ||S.size1()!=S.size2() || J.size1()!=S.size1() ){
        cerr << " \n Loewdin transformation only works for quadratic matrices " << endl;
        return -1;
    }
  
    ub::vector<double> S_eigenvalues;
    linalg_eigenvalues( S_eigenvalues, S);
    if ( S_eigenvalues[0] < 0.0 ) {
        cerr << " \n Negative eigenvalues in Loewdin transformation " << endl;
        return -1;
    }

    ub::matrix<double> _diagS = ub::zero_matrix<double>(J.size1(),J.size2() );
     for ( unsigned _i =0; _i < J.size1() ; _i++){

         _diagS(_i,_i) = 1.0/sqrt(S_eigenvalues[_i]);
     }
    ub::matrix<double> _temp = ub::prod( _diagS, ub::trans(S));
    S = ub::prod( S,_temp );
    
    _temp = ub::prod(J, S);
    J = ub::prod( S, _temp);
    return S_eigenvalues[0];
    
}



int linalg_matrixsqrt(ub::matrix<double> &S){
    if (S.size1()!=S.size2()) {
        cerr << " \nMatrix sqrt only works for quadratic matrices " << endl;
        return -1;
    }
  
    ub::vector<double> S_eigenvalues;
    linalg_eigenvalues( S_eigenvalues, S);
    if ( S_eigenvalues[0] < 0.0 ) {
        cerr << " \n Negative eigenvalues in matrix_sqrt transformation " << endl;
        return -1;
    }

    ub::matrix<double> _diagS = ub::zero_matrix<double>(S.size1(),S.size2() );
     for ( unsigned _i =0; _i < S.size1() ; _i++){

         _diagS(_i,_i) = sqrt(S_eigenvalues[_i]);
     }
    ub::matrix<double> _temp = ub::prod( _diagS, ub::trans(S));
    S = ub::prod( S,_temp );
    
  
    return 1;
} 
    
    double linalg_getMax( const ub::matrix<double>& _matrix, bool absolut ){

   double _maximum = 0.0;
   if(absolut){
   for ( unsigned _i = 0; _i < _matrix.size1(); _i++ ){
     for ( unsigned _j = 0; _j < _matrix.size2(); _j++ ){
       if ( std::abs(_matrix(_i,_j)) > _maximum ) {
	 _maximum = std::abs(_matrix(_i,_j));
       }				   
     }
   }
   }
   else{
      for ( unsigned _i = 0; _i < _matrix.size1(); _i++ ){
     for ( unsigned _j = 0; _j < _matrix.size2(); _j++ ){
       if ( _matrix(_i,_j) > _maximum ) {
	 _maximum = _matrix(_i,_j);
       }				   
     }
   } 
   }
   return _maximum;
   }
    
    
    ub::matrix<double> linalg_abs( const ub::matrix<double>& _matrix ){

        ub::matrix<double> A=ub::matrix<double>(_matrix.size1(),_matrix.size2());

   for ( unsigned _i = 0; _i < _matrix.size1(); _i++ ){
     for ( unsigned _j = 0; _j < _matrix.size2(); _j++ ){
         A(_i,_j)=std::abs(_matrix(_i,_j));
       }				   
     }
   

   return A;
   }
    
   double linalg_getRMS(const ub::matrix<double>& _matrix ){

   double _rms = 0.0;
   int _n = 0;

   for ( unsigned _i = 0; _i < _matrix.size1(); _i++ ){
     for ( unsigned _j = 0; _j < _matrix.size2(); _j++ ){
       _rms += _matrix(_i,_j) * _matrix(_i,_j);
       _n++;
     }
   }

   _rms = sqrt(_rms/_n);

   return _rms;
   }
   
   double linalg_traceofProd(const ub::matrix<double>& A,const ub::matrix<double>& B ){
       double trace=0.0;
       const ub::vector<double>& Aarray=A.data();
       const ub::vector<double>& Barray=B.data();
        if (Aarray.size()!=Barray.size() || A.size1()!=B.size2() || A.size2()!=B.size1()) {
        cerr << " \ntraceofProd Matrcies have different number of indexes " << endl;
        return -1;
       }
       for (unsigned i=0;i<Aarray.size();i++){
           trace+=Aarray(i)*Barray(i);
       }
       return trace;
   }
   
   
   ub::vector<double> linalg_abs(const ub::vector<double>& a){
       ub::vector<double> b=ub::zero_vector<double>(a.size());
       for (unsigned i=0;i<a.size();i++){
           b(i)=std::abs(a(i));
       }
               
       return b;
   }
   
    double linalg_getMax( const ub::vector<double>& a, bool absolut ){

   double _maximum = 0.0;
   if(absolut){
   for ( unsigned _i = 0; _i < a.size(); _i++ ){
       if ( std::abs(a(_i)) > _maximum ) {
	 _maximum = std::abs(a(_i));
       }				   
     
   }
   }
   else{
     for ( unsigned _i = 0; _i < a.size(); _i++ ){
       if ( a(_i) > _maximum ) {
	 _maximum = a(_i);
       }				   
     
   }
   }
   return _maximum;
   }
 
    



}}

