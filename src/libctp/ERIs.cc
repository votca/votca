/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICEN_olE-2.0
 *
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/ctp/aomatrix.h>
#include <votca/ctp/aobasis.h>
#include <votca/ctp/threecenters.h>
#include <votca/ctp/ERIs.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <votca/tools/linalg.h>
#include <omp.h>

using namespace std;
using namespace votca::tools;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;
        
        
        
        void ERIs::Initialize (AOBasis &_dftbasis, AOBasis &_auxbasis,  AOOverlap &_auxAOverlap, AOCoulomb &_auxAOcoulomb){

           
           
            _threecenter.Fill( _auxbasis, _dftbasis );
            
            ub::matrix<double> _inverse=ub::zero_matrix<double>( _auxAOverlap.Dimension(), _auxAOverlap.Dimension());
            //_auxAOverlap.Print("auxAO");
            linalg_invert( _auxAOverlap.Matrix() , _inverse);
            /*
            for (int j=0;j<_inverse.size2();j++){
             for (int i=0;i<_inverse.size1();i++){
                cout << "_inverse ("<< i <<":"<< j<<")="<<_inverse(i,j)<<endl;
             }}
            exit(0);
            */
            ub::matrix<double> _test=ub::prod(_auxAOverlap.Matrix(),_inverse);
            /*
            for (int i=0;i<_test.size1();i++){
                cout << "_test ("<< i <<")="<<_test(i,i)<<endl;
            }
            exit(0);
            */
            ub::matrix<double> _temp=ub::prod(_auxAOcoulomb.Matrix(),_inverse);
            _Vcoulomb=ub::prod(_inverse,_temp);
            /*
              for (int j=0;j<_Vcoulomb.size2();j++){
             for (int i=0;i<_Vcoulomb.size1();i++){
                cout << "_Vcoulomb ("<< i <<":"<< j<<")="<<_Vcoulomb(i,j)<<endl;
             }}
            exit(0);
             */
            //cout << endl;
            //cout << _Vcoulomb.size1() << "x" << _Vcoulomb.size2()<<"Vsize"<< endl;
            //cout << _auxbasis.AOBasisSize() << " aobasissize"<<endl;

        
        
        }
        
        
        void ERIs::CalculateERIs (ub::matrix<double> &DMAT){
            _ERIs=ub::zero_matrix<double>(DMAT.size1(),DMAT.size2());
            ub::vector<double> dmatasarray=DMAT.data();
            ub::matrix<double> Itilde=ub::zero_matrix<double>(_threecenter.getSize(),1);
            //cout << _threecenter.getSize() << " Size-Threecenter"<<endl;
            //check Efficiency !!!! someday 
            for ( int _i=0; _i<_threecenter.getSize();_i++){
                ub::vector<double>threecenterasarray=(_threecenter.getDatamatrix(_i)).data();
                //cout << _threecenter.getDatamatrix(_i).size1() << "x"<< _threecenter.getDatamatrix(_i).size2() <<" Size-Threecenter,matrix"<<endl;
                // Trace over prod::DMAT,I(l)=componentwise product over 
                for ( int _j=0; _j<threecenterasarray.size();_j++){
                    Itilde(_i,0)+=dmatasarray[_j]*threecenterasarray[_j];
                }
            }
            //cout << "Itilde " <<Itilde << endl;
            ub::matrix<double>K=ub::prod(_Vcoulomb,Itilde);
            //cout << "K " << K << endl;
            for ( int _i=0; _i<K.size1(); _i++){
                
            _ERIs+=_threecenter.getDatamatrix(_i)*K(_i,0);    
            //cout << "I " << _threecenter.getDatamatrix(_i) << endl;
            //cout<< "ERIs " <<_ERIs<< endl;
            }
            
           
            CalculateEnergy(dmatasarray);
        }
        
        
        
        
        
        // This is simply Trace(prod(ERIs*Dmat))=Totenergy
        void ERIs::CalculateEnergy(ub::vector<double> &dmatasarray){
            _ERIsenergy=0;
            ub::vector<double> ERIsasarray=_ERIs.data();
            for ( int _i=0;_i<ERIsasarray.size();_i++){
                _ERIsenergy+=dmatasarray[_i]*ERIsasarray[_i];
                
            }
            
            
            
            
        }
        
        
        void ERIs::printERIs(){
          for (int i=0; i< _ERIs.size1(); i++){
                for (int j=0; j< _ERIs.size2();j++){
                    cout << "ERIs [" << i<<":"<<j<<"]="<<_ERIs(i,j)<<endl;
                }
            }
        }
        
        
        
        
        
    }}