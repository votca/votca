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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/ERIs.h>
#include <votca/xtp/symmetric_matrix.h>

namespace votca {
    namespace xtp {
        
        
        
        
    void ERIs::Initialize(AOBasis &_dftbasis, AOBasis &_auxbasis,const Eigen::MatrixXd &inversesqrt_Coulomb) {
            _threecenter.Fill( _auxbasis, _dftbasis,inversesqrt_Coulomb);
            return;
        }

        void ERIs::Initialize_4c_small_molecule(AOBasis &_dftbasis) {
          _fourcenter.Fill_4c_small_molecule( _dftbasis );
          return;
        }
        
        
void ERIs::CalculateERIs(const Eigen::MatrixXd &DMAT) {
      

      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<Eigen::MatrixXd >ERIS_thread;

      for (unsigned i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd thread = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
        ERIS_thread.push_back(thread);
      }

#pragma omp parallel for
      for (unsigned thread = 0; thread < nthreads; ++thread) {
        Symmetric_Matrix dmat_sym = Symmetric_Matrix(DMAT);
        for (int _i = thread; _i < _threecenter.getSize(); _i += nthreads) {
          const Symmetric_Matrix &threecenter = _threecenter.getDatamatrix(_i);
          // Trace over prod::DMAT,I(l)=componentwise product over 
          double factor = threecenter.TraceofProd(dmat_sym);
          threecenter.AddtoEigenMatrix(ERIS_thread[thread], factor);
        }
      }
      
      _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
      for (const auto& thread : ERIS_thread) {
        _ERIs += thread;
      }

      CalculateEnergy(DMAT);
      return;
    }

    void ERIs::CalculateEXX(const Eigen::MatrixXd &DMAT) {
      
      int nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<Eigen::MatrixXd >EXX_thread;

      for (int i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd thread = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
        EXX_thread.push_back(thread);
      }
      
      #pragma omp parallel for
      for (int thread = 0; thread < nthreads; ++thread) {
        Eigen::MatrixXd D=DMAT;
        for(int i=thread;i<_threecenter.getSize();i+= nthreads){
          const Eigen::MatrixXd threecenter = _threecenter.getDatamatrix(i).FullMatrix();
          EXX_thread[thread]+=threecenter*D*threecenter;
        }
      }
      _EXXs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
      for (const auto& thread : EXX_thread) {
        _EXXs += thread;
      }

      CalculateEXXEnergy(DMAT);
      return;
    }
    
     void ERIs::CalculateEXX(const Eigen::Block<Eigen::MatrixXd>& occMos,const Eigen::MatrixXd  &DMAT) {
      
      int nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<Eigen::MatrixXd >EXX_thread;

      for (int i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd thread = Eigen::MatrixXd::Zero(occMos.rows(), occMos.rows());
        EXX_thread.push_back(thread);
      }
      
      #pragma omp parallel for
      for (int thread = 0; thread < nthreads; ++thread) {
        Eigen::MatrixXd occ=occMos;
        for(int i=thread;i<_threecenter.getSize();i+= nthreads){
          const Eigen::MatrixXd TCxMOs_T = occ.transpose()*_threecenter.getDatamatrix(i).FullMatrix();
          EXX_thread[thread]+=TCxMOs_T.transpose()*TCxMOs_T;
        }
      }
      _EXXs = Eigen::MatrixXd::Zero(occMos.rows(), occMos.rows());
      for (const auto& thread : EXX_thread) {
        _EXXs += 2*thread;
      }

      CalculateEXXEnergy(DMAT);
      return;
    }
    



        void ERIs::CalculateERIs_4c_small_molecule(const Eigen::MatrixXd  &DMAT) {

          _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
          
          const Eigen::VectorXd& _4c_vector = _fourcenter.get_4c_vector();
          
          

          int dftBasisSize = DMAT.rows();
          int vectorSize = (dftBasisSize*(dftBasisSize+1))/2;
          #pragma omp parallel for
          for (int _i = 0; _i < dftBasisSize; _i++) {
            int sum_i = (_i*(_i+1))/2;
            for (int _j = _i; _j < dftBasisSize; _j++) {
              int _index_ij = dftBasisSize * _i - sum_i + _j;
              int _index_ij_kl_a = vectorSize * _index_ij - (_index_ij*(_index_ij+1))/2;
              for (int _k = 0; _k < dftBasisSize; _k++) {
                int sum_k = (_k*(_k+1))/2;
                for (int _l = _k; _l < dftBasisSize; _l++) {
                  int _index_kl = dftBasisSize * _k - sum_k + _l;

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

          CalculateEnergy(DMAT);
          return;
        }
        
        
        void ERIs::CalculateEXX_4c_small_molecule(const Eigen::MatrixXd &DMAT) {
            int nthreads = 1;
            #ifdef _OPENMP
                  nthreads = omp_get_max_threads();
            #endif  
          
          
          std::vector<Eigen::MatrixXd >EXX_thread;

      for (int i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd thread = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
        EXX_thread.push_back(thread);
      }
          
          const Eigen::VectorXd& _4c_vector = _fourcenter.get_4c_vector();

          int dftBasisSize = DMAT.rows();
          int vectorSize = (dftBasisSize*(dftBasisSize+1))/2;
        #pragma omp parallel for
        for (int thread = 0; thread < nthreads; ++thread) {
          for (int _i = thread; _i < dftBasisSize; _i+= nthreads) {
            int sum_i = (_i*(_i+1))/2;
            for (int _j = _i; _j < dftBasisSize; _j++) {
              int _index_ij = DMAT.cols() * _i - sum_i + _j;
              int _index_ij_kl_a = vectorSize * _index_ij - (_index_ij*(_index_ij+1))/2;
              for (int _k = 0; _k < dftBasisSize; _k++) {
                int sum_k = (_k*(_k+1))/2;
                for (int _l = _k; _l < dftBasisSize; _l++) {
                  int _index_kl = DMAT.cols() * _k - sum_k + _l;

                  int _index_ij_kl = _index_ij_kl_a + _index_kl;
                  if (_index_ij > _index_kl) _index_ij_kl = vectorSize * _index_kl - (_index_kl*(_index_kl+1))/2 + _index_ij;
                  double factorij=1;
                  if(_i==_j){factorij=0.5;}
                  double factorkl=1;
                  if(_l==_k){factorkl=0.5;}
                  double factor=factorij*factorkl;
                  EXX_thread[thread](_i, _l) +=factor*DMAT(_j, _k) * _4c_vector(_index_ij_kl);
                  EXX_thread[thread](_j, _l) +=factor*DMAT(_i, _k) * _4c_vector(_index_ij_kl);
                  EXX_thread[thread](_i, _k) +=factor*DMAT(_j, _l) * _4c_vector(_index_ij_kl);
                  EXX_thread[thread](_j, _k) +=factor*DMAT(_i, _l) * _4c_vector(_index_ij_kl);
                }
              }
            }
          }
        }
          _EXXs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
          for (const auto& thread : EXX_thread) {
          _EXXs += thread;
          }

          CalculateEXXEnergy(DMAT);
          return;
        }
        
        
        
        
         void ERIs::CalculateEnergy(const Eigen::MatrixXd &DMAT){
            _ERIsenergy=_ERIs.cwiseProduct(DMAT).sum();
            return;
        }
         
         
         void ERIs::CalculateEXXEnergy(const Eigen::MatrixXd &DMAT){
           _EXXsenergy=_EXXs.cwiseProduct(DMAT).sum();
            return;
        }
        
        
        void ERIs::printERIs(){
          for (int i=0; i< _ERIs.cols(); i++){
                for (int j=0; j< _ERIs.rows();j++){
                    cout << "ERIs [" << i<<":"<<j<<"]="<<_ERIs(i,j)<<endl;
                }
            }
          return;
        }
        
        
        
        
        
    }}
