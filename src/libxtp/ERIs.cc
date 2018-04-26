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
        
        
        
        
    void ERIs::Initialize(AOBasis &_dftbasis, AOBasis &_auxbasis,const Eigen::MatrixXd &inverse_Coulomb) {

          _inverse_Coulomb=inverse_Coulomb;
          _threecenter.Fill( _auxbasis, _dftbasis );
          return;
        }



        void ERIs::Initialize_4c_small_molecule(AOBasis &_dftbasis) {
          _fourcenter.Fill_4c_small_molecule( _dftbasis );
          return;
        }

        
        
        void ERIs::Initialize_4c_screening(AOBasis &_dftbasis, double eps) {
          
          _with_screening = true;
          _screening_eps = eps;
          CalculateERIsDiagonals(_dftbasis);
          return;
        }
        
        
        
        void ERIs::CalculateERIs (const Eigen::MatrixXd &DMAT){
          
           Symmetric_Matrix dmat_sym=Symmetric_Matrix(DMAT);
            _ERIs=Eigen::MatrixXd::Zero(DMAT.rows(),DMAT.cols());
            Eigen::VectorXd Itilde=Eigen::VectorXd::Zero(_threecenter.getSize());
          
            #pragma omp parallel for
            for ( int _i=0; _i<_threecenter.getSize();_i++){
                const Symmetric_Matrix &threecenter=_threecenter.getDatamatrix(_i);
                // Trace over prod::DMAT,I(l)=componentwise product over 
                
                Itilde(_i)=threecenter.TraceofProd(dmat_sym);
            }
            const Eigen::VectorXd K=_inverse_Coulomb*Itilde;
            
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               std::vector<Eigen::MatrixXd >ERIS_thread;
               
               for(unsigned i=0;i<nthreads;++i){
                   Eigen::MatrixXd thread=Eigen::MatrixXd::Zero(_ERIs.rows(),_ERIs.cols());
                   ERIS_thread.push_back(thread);
               }
            
            #pragma omp parallel for
            for (unsigned thread=0;thread<nthreads;++thread){
                for ( unsigned _i = thread; _i < K.size(); _i+=nthreads){
                _threecenter.getDatamatrix(_i).AddtoEigenMatrix(ERIS_thread[thread],K(_i));    
                }
            }
              
            for (unsigned thread=0;thread<nthreads;++thread){
                _ERIs+=ERIS_thread[thread];
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
        
        
        void ERIs::CalculateERIs_4c_direct(const AOBasis& dftbasis, const Eigen::MatrixXd &DMAT) {

          // Number of shells
          int numShells = dftbasis.getNumofShells();
          
          // Initialize ERIs matrix
          _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

          #pragma omp parallel
          { // Begin omp parallel
            
            
            Eigen::MatrixXd ERIs_thread = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
            
            #pragma omp for
            for (int iShell_3 = 0; iShell_3 < numShells; iShell_3++) {
              const AOShell& shell_3 = *dftbasis.getShell(iShell_3);
              for (int iShell_4 = iShell_3; iShell_4 < numShells; iShell_4++) {
                const AOShell& shell_4 = *dftbasis.getShell(iShell_4);
                for (int iShell_1 = iShell_3; iShell_1 < numShells; iShell_1++) {
                  const AOShell& shell_1 = *dftbasis.getShell(iShell_1);
                  for (int iShell_2 = iShell_1; iShell_2 < numShells; iShell_2++) {
                    const AOShell& shell_2 = *dftbasis.getShell(iShell_2);

                    // Pre-screening
                    if (_with_screening && CheckScreen(_screening_eps, shell_1, shell_2, shell_3, shell_4))
                      continue;

                    // Get the current 4c block
                    Eigen::MatrixXd subMatrix = Eigen::MatrixXd::Zero(shell_1.getNumFunc() * shell_2.getNumFunc(), shell_3.getNumFunc() * shell_4.getNumFunc());
                    bool nonzero = _fourcenter.FillFourCenterRepBlock(subMatrix, &shell_1, &shell_2, &shell_3, &shell_4);

                    // If there are only zeros, we don't need to put anything in the ERIs matrix
                    if (!nonzero)
                      continue;

                    // Begin fill ERIs matrix

                    FillERIsBlock(ERIs_thread, DMAT, subMatrix, shell_1, shell_2, shell_3, shell_4);

                    // Symmetry 1 <--> 2
                    if (iShell_1 != iShell_2)
                      FillERIsBlock(ERIs_thread, DMAT, subMatrix, shell_2, shell_1, shell_3, shell_4);

                    // Symmetry 3 <--> 4
                    if (iShell_3 != iShell_4)
                      FillERIsBlock(ERIs_thread, DMAT, subMatrix, shell_1, shell_2, shell_4, shell_3);

                    // Symmetry 1 <--> 2 and 3 <--> 4
                    if (iShell_1 != iShell_2 && iShell_3 != iShell_4)
                      FillERIsBlock(ERIs_thread, DMAT, subMatrix, shell_2, shell_1, shell_4, shell_3);

                    // Symmetry (1, 2) <--> (3, 4)
                    if (iShell_1 != iShell_3) {

                      // We need the transpose of "subMatrix"
                      Eigen::MatrixXd subMatrix2 = subMatrix.transpose();

                      FillERIsBlock(ERIs_thread, DMAT, subMatrix2, shell_3, shell_4, shell_1, shell_2);

                      // Symmetry 1 <--> 2
                      if (iShell_1 != iShell_2)
                        FillERIsBlock(ERIs_thread, DMAT, subMatrix2, shell_3, shell_4, shell_2, shell_1);

                      // Symmetry 3 <--> 4
                      if (iShell_3 != iShell_4)
                        FillERIsBlock(ERIs_thread, DMAT, subMatrix2, shell_4, shell_3, shell_1, shell_2);

                      // Symmetry 1 <--> 2 and 3 <--> 4
                      if (iShell_1 != iShell_2 && iShell_3 != iShell_4)
                        FillERIsBlock(ERIs_thread, DMAT, subMatrix2, shell_4, shell_3, shell_2, shell_1);
                    }

                    // End fill ERIs matrix
                  } // End loop over shell 2
                } // End loop over shell 1
              } // End loop over shell 4
            } // End loop over shell 3
            
            #pragma omp critical
            {    
              _ERIs += ERIs_thread;
            }
          } 

          // Fill lower triangular part using symmetry
          for (size_t i = 0; i < DMAT.cols(); i++){
            for (size_t j = i + 1; j < DMAT.rows(); j++){
              _ERIs(j, i) = _ERIs(i, j);
            }
          }

          CalculateEnergy(DMAT);
          return;
        }
        
        
        void ERIs::FillERIsBlock(Eigen::MatrixXd& ERIsCur, const Eigen::MatrixXd& DMAT,
                const Eigen::MatrixXd& subMatrix,
                const AOShell& shell_1, const AOShell& shell_2,
                const AOShell& shell_3, const AOShell& shell_4) {

          for (int iFunc_3 = 0; iFunc_3 < shell_3.getNumFunc(); iFunc_3++) {
            int ind_3 = shell_3.getStartIndex() + iFunc_3;
            for (int iFunc_4 = 0; iFunc_4 < shell_4.getNumFunc(); iFunc_4++) {
              int ind_4 = shell_4.getStartIndex() + iFunc_4;

              // Symmetry
              if (ind_3 > ind_4)
                continue;

              // Column index in the current sub-matrix
              int ind_subm_34 = shell_3.getNumFunc() * iFunc_4 + iFunc_3;

              for (int iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
                int ind_1 = shell_1.getStartIndex() + iFunc_1;
                for (int iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
                  int ind_2 = shell_2.getStartIndex() + iFunc_2;

                  // Symmetry
                  if (ind_1 > ind_2)
                    continue;

                  // Row index in the current sub-matrix
                  int ind_subm_12 = shell_1.getNumFunc() * iFunc_2 + iFunc_1;

                  //Symmetry for diagonal elements
                  double multiplier=(ind_1 == ind_2 ? 1.0 : 2.0);
                  // Fill ERIs matrix
                  ERIsCur(ind_3, ind_4) += multiplier* DMAT(ind_1, ind_2) * subMatrix(ind_subm_12, ind_subm_34);
                } // End loop over functions in shell 2
              } // End loop over functions in shell 1
            } // End loop over functions in shell 4
          } // End loop over functions in shell 3
          
          return;
        }
        
        
        void ERIs::CalculateERIsDiagonals(const AOBasis& dftbasis) {
          
          // Number of shells
          int numShells = dftbasis.getNumofShells();
          // Total number of functions
          int dftBasisSize = dftbasis.AOBasisSize();
          
          _diagonals = Eigen::MatrixXd::Zero(dftBasisSize,dftBasisSize);
          
          for (int iShell_1 = 0; iShell_1 < numShells; iShell_1++) {
            const AOShell& shell_1 = *dftbasis.getShell(iShell_1);
            for (int iShell_2 = iShell_1; iShell_2 < numShells; iShell_2++) {
              const AOShell& shell_2 = *dftbasis.getShell(iShell_2);
              
              // Get the current 4c block
              Eigen::MatrixXd subMatrix = Eigen::MatrixXd::Zero(shell_1.getNumFunc() * shell_2.getNumFunc(), shell_1.getNumFunc() * shell_2.getNumFunc());
              bool nonzero = _fourcenter.FillFourCenterRepBlock(subMatrix, &shell_1, &shell_2, &shell_1, &shell_2);
              
              if (!nonzero)
                continue;
              
              int index = 0; // Diagonal index
              
              for (int iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
                int ind_1 = shell_1.getStartIndex() + iFunc_1;
                for (int iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
                  int ind_2 = shell_2.getStartIndex() + iFunc_2;

                  // Symmetry
                  if (ind_1 > ind_2)
                    continue;
    
                  _diagonals(ind_1, ind_2) = subMatrix(index, index);
                  
                  // Symmetry
                  if (ind_1 != ind_2){
                    _diagonals(ind_2, ind_1) = _diagonals(ind_1, ind_2);
                  }
                 
                  index++; // composite index of shell1 and shell2
                }
              }
            }
          }
          
          return;
        }
        

        bool ERIs::CheckScreen(double eps, const AOShell& shell_1, const AOShell& shell_2, const AOShell& shell_3, const AOShell& shell_4) {
          const double eps2=eps * eps;
          
          for (int iFunc_3 = 0; iFunc_3 < shell_3.getNumFunc(); iFunc_3++) {
            int ind_3 = shell_3.getStartIndex() + iFunc_3;
            for (int iFunc_4 = 0; iFunc_4 < shell_4.getNumFunc(); iFunc_4++) {
              int ind_4 = shell_4.getStartIndex() + iFunc_4;

              // Symmetry
              if (ind_3 > ind_4){
                continue;
              }
              
              for (int iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
                int ind_1 = shell_1.getStartIndex() + iFunc_1;
                for (int iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
                  int ind_2 = shell_2.getStartIndex() + iFunc_2;

                  // Symmetry
                  if (ind_1 > ind_2){
                    continue;
                  }
                  // Cauchy–Schwarz
                  // <ab|cd> <= sqrt(<ab|ab>) * sqrt(<cd|cd>)
                  double ub = _diagonals(ind_1, ind_2) * _diagonals(ind_3, ind_4);
                  
                  // Compare with tolerance
                  if (ub > eps2){
                    return false; // We must compute ERIS for the whole block
                  }
                }
              }
            }
          }
          
          return true; // We can skip the whole block
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
