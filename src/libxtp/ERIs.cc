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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/ERIs.h>



using namespace votca::tools;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
        
        
        
    void ERIs::Initialize(AOBasis &_dftbasis, AOBasis &_auxbasis) {

           
           
            _threecenter.Fill( _auxbasis, _dftbasis );
            //cout << "_threeceenter(0)"<<endl;
            //cout << _threecenter.getDatamatrix(0)<<endl;

            return;
        }



        void ERIs::Initialize_4c_small_molecule(AOBasis &_dftbasis) {

          _fourcenter.Fill_4c_small_molecule( _dftbasis );

          return;
        }

         
        
        

  /*  void Eris::Initialize_Symmetric(AOBasis &_dftbasis,AOBasis &_auxbasis, AOCoulomb &_auxcoulomb){
        _threecenter.Fill(_auxbasis,_dftbasis);
        eigenlinalg_matrixsqrt(_auxcoulomb.Matrix());
        
        for 
        
    }    
        
 */       
        
        
      
        
        
        void ERIs::CalculateERIs (const ub::matrix<double> &DMAT, const ub::matrix<double> &_inverse_Coulomb){

            //cout << _auxAOcoulomb.Matrix()<<endl;
            //cout << "inverse Coulomb"<< endl;
            //cout << _inverse_Coulomb<<endl;
            
            _ERIs=ub::zero_matrix<double>(DMAT.size1());
           
        
            const ub::vector<double>& dmatasarray=DMAT.data();
          
            ub::matrix<double> Itilde=ub::zero_matrix<double>(_threecenter.getSize(),1);
            //cout << _threecenter.getSize() << " Size-Threecenter"<<endl;
            //check Efficiency !!!! someday 
            #pragma omp parallel for
            for ( int _i=0; _i<_threecenter.getSize();_i++){
                ub::symmetric_matrix<double> &threecenter=_threecenter.getDatamatrix(_i);
                // Trace over prod::DMAT,I(l)=componentwise product over 
                for ( unsigned _j=0; _j<DMAT.size1();_j++){
                    for(unsigned _k=0;_k<DMAT.size2();_k++){
                    Itilde(_i,0)+=DMAT(_j,_k)*threecenter(_j,_k);
                    }
                }
            }
            //cout << "Itilde " <<Itilde << endl;
            ub::matrix<double>K=ub::prod(_inverse_Coulomb,Itilde);
            //cout << "K " << K << endl;
            for ( unsigned _i = 0; _i < K.size1(); _i++){
                
            _ERIs+=_threecenter.getDatamatrix(_i)*K(_i,0);    
            //cout << "I " << _threecenter.getDatamatrix(_i) << endl;
            //cout<< "ERIs " <<_ERIs<< endl;
            }

            CalculateEnergy(dmatasarray);
            //cout<<_ERIs<<endl;
        }














        void ERIs::CalculateERIs_4c_small_molecule(const ub::matrix<double> &DMAT) {

          _ERIs = ub::zero_matrix<double>(DMAT.size1(), DMAT.size2());
          const ub::vector<double> dmatasarray = DMAT.data();
          const ub::vector<double>& _4c_vector = _fourcenter.get_4c_vector();

          int dftBasisSize = DMAT.size1();
          int vectorSize = (dftBasisSize*(dftBasisSize+1))/2;
          #pragma omp parallel for
          for (unsigned _i = 0; _i < DMAT.size1(); _i++) {
            unsigned sum_i = (_i*(_i+1))/2;
            for (unsigned _j = _i; _j < DMAT.size2(); _j++) {
              unsigned _index_ij = DMAT.size2() * _i - sum_i + _j;
              unsigned _index_ij_kl_a = vectorSize * _index_ij - (_index_ij*(_index_ij+1))/2;
              for (unsigned _k = 0; _k < DMAT.size1(); _k++) {
                unsigned sum_k = (_k*(_k+1))/2;
                for (unsigned _l = _k; _l < DMAT.size2(); _l++) {
                  unsigned _index_kl = DMAT.size2() * _k - sum_k + _l;

                  unsigned _index_ij_kl = _index_ij_kl_a + _index_kl;
                  if (_index_ij > _index_kl) _index_ij_kl = vectorSize * _index_kl - (_index_kl*(_index_kl+1))/2 + _index_ij;

                  if (_l == _k) {
                    _ERIs(_i, _j) += DMAT(_k, _l) * _4c_vector(_index_ij_kl);
                  } else {
                    _ERIs(_i, _j) += 2. * DMAT(_k, _l) * _4c_vector(_index_ij_kl);
                  }

                }
              }
              _ERIs(_j, _i) = _ERIs(_i, _j);
            }
          }

          CalculateEnergy(dmatasarray);

        }
        
    








        
         void ERIs::CalculateEnergy(const ub::vector<double> &dmatasarray){
            
            const ub::vector<double>& ERIsasarray=_ERIs.data();
            double energy=0.0;
           #pragma omp parallel for reduction(+:energy)
            for (unsigned _i=0;_i<ERIsasarray.size();_i++){
                energy+=dmatasarray[_i]*ERIsasarray[_i];
                
            }
            _ERIsenergy=energy;

        }
        
        
        void ERIs::printERIs(){
          for (unsigned i=0; i< _ERIs.size1(); i++){
                for (unsigned j=0; j< _ERIs.size2();j++){
                    cout << "ERIs [" << i<<":"<<j<<"]="<<_ERIs(i,j)<<endl;
                }
            }
        }
        
        
        
        
        
    }}
