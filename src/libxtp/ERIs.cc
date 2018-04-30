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
      Symmetric_Matrix dmat_sym = Symmetric_Matrix(DMAT);
      _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<Eigen::MatrixXd >ERIS_thread;

      for (unsigned i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd thread = Eigen::MatrixXd::Zero(_ERIs.rows(), _ERIs.cols());
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

      for (const auto& thread : ERIS_thread) {
        _ERIs += thread;
      }

      CalculateEnergy(DMAT);
      return;
    }




        void ERIs::CalculateERIs_4c_small_molecule(const Eigen::MatrixXd  &DMAT) {

          _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
          
          const Eigen::VectorXd& _4c_vector = _fourcenter.get_4c_vector();

          int dftBasisSize = DMAT.rows();
          int vectorSize = (dftBasisSize*(dftBasisSize+1))/2;
          #pragma omp parallel for
          for (unsigned _i = 0; _i < DMAT.rows(); _i++) {
            unsigned sum_i = (_i*(_i+1))/2;
            for (unsigned _j = _i; _j < DMAT.cols(); _j++) {
              unsigned _index_ij = DMAT.cols() * _i - sum_i + _j;
              unsigned _index_ij_kl_a = vectorSize * _index_ij - (_index_ij*(_index_ij+1))/2;
              for (unsigned _k = 0; _k < DMAT.rows(); _k++) {
                unsigned sum_k = (_k*(_k+1))/2;
                for (unsigned _l = _k; _l < DMAT.cols(); _l++) {
                  unsigned _index_kl = DMAT.cols() * _k - sum_k + _l;

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

          _EXXs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
          
          const Eigen::VectorXd& _4c_vector = _fourcenter.get_4c_vector();

          int dftBasisSize = DMAT.rows();
          int vectorSize = (dftBasisSize*(dftBasisSize+1))/2;
          #pragma omp parallel for
          for (unsigned _i = 0; _i < DMAT.rows(); _i++) {
            unsigned sum_i = (_i*(_i+1))/2;
            for (unsigned _j = _i; _j < DMAT.cols(); _j++) {
              unsigned _index_ij = DMAT.cols() * _i - sum_i + _j;
              unsigned _index_ij_kl_a = vectorSize * _index_ij - (_index_ij*(_index_ij+1))/2;
              for (unsigned _k = 0; _k < DMAT.rows(); _k++) {
                unsigned sum_k = (_k*(_k+1))/2;
                for (unsigned _l = _k; _l < DMAT.cols(); _l++) {
                  unsigned _index_kl = DMAT.cols() * _k - sum_k + _l;

                  unsigned _index_ij_kl = _index_ij_kl_a + _index_kl;
                  if (_index_ij > _index_kl) _index_ij_kl = vectorSize * _index_kl - (_index_kl*(_index_kl+1))/2 + _index_ij;

                  if (_l == _k) {
                    _EXXs(_i, _l) += DMAT(_j, _k) * _4c_vector(_index_ij_kl);
                  } else {
                    _EXXs(_i, _l) += 2. * DMAT(_j, _k) * _4c_vector(_index_ij_kl);
                  }

                }
              }
              _EXXs(_j, _i) = _EXXs(_i, _j);
            }
          }

          CalculateEXXEnergy(DMAT);
          return;
        }
        
        
        
        
         void ERIs::CalculateEnergy(const Eigen::MatrixXd &DMAT){
            _ERIsenergy=_ERIs.cwiseProduct(DMAT).sum();
            return;
        }
         
         
         void ERIs::CalculateEXXEnergy(const Eigen::MatrixXd &DMAT){
           _ERIsenergy=_EXXs.cwiseProduct(DMAT).sum();
            return;
        }
        
        
        void ERIs::printERIs(){
          for (unsigned i=0; i< _ERIs.cols(); i++){
                for (unsigned j=0; j< _ERIs.rows();j++){
                    cout << "ERIs [" << i<<":"<<j<<"]="<<_ERIs(i,j)<<endl;
                }
            }
          return;
        }
        
        
        
        
        
    }}
