/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

        void ERIs::Initialize(AOBasis &dftbasis, AOBasis &auxbasis,const Eigen::MatrixXd &inversesqrt_Coulomb) {
          
          _threecenter.Fill( auxbasis, dftbasis,inversesqrt_Coulomb);
          return;
        }

        
        void ERIs::Initialize_4c_small_molecule(AOBasis &dftbasis) {
          
          _fourcenter.Fill_4c_small_molecule( dftbasis );
          return;
        }

        
        void ERIs::Initialize_4c_screening(AOBasis &dftbasis, double eps) {
          
          _with_screening = true;
          _screening_eps = eps;
          CalculateERIsDiagonals(dftbasis);
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
            for (int i = thread; i < _threecenter.size(); i += nthreads) {
              const Symmetric_Matrix &threecenter = _threecenter.getDatamatrix(i);
              // Trace over prod::DMAT,I(l)=componentwise product over 
              const double factor = threecenter.TraceofProd(dmat_sym);
              Eigen::SelfAdjointView<Eigen::MatrixXd,Eigen::Upper> m=ERIS_thread[thread].selfadjointView<Eigen::Upper>();
              threecenter.AddtoEigenUpperMatrix(m, factor);
            }
          }

          _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
          for (const auto& thread : ERIS_thread) {
            _ERIs += thread;
          }
          _ERIs+=_ERIs.triangularView<Eigen::StrictlyUpper>().transpose();

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
            for(int i=thread;i<_threecenter.size();i+= nthreads){
              const Eigen::MatrixXd threecenter = _threecenter.getDatamatrix(i).UpperMatrix();
              EXX_thread[thread].noalias()+=threecenter.selfadjointView<Eigen::Upper>()*D*threecenter.selfadjointView<Eigen::Upper>();
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
            for(int i=thread;i<_threecenter.size();i+= nthreads){
              const Eigen::MatrixXd TCxMOs_T = occ.transpose()*_threecenter.getDatamatrix(i).UpperMatrix().selfadjointView<Eigen::Upper>();
              EXX_thread[thread].noalias()+=TCxMOs_T.transpose()*TCxMOs_T;
            }
          }
          _EXXs = Eigen::MatrixXd::Zero(occMos.rows(), occMos.rows());
          for (const auto& thread : EXX_thread) {
            _EXXs += 2*thread;
          }

          CalculateEXXEnergy(DMAT);
          return;
        }

        
        void ERIs::CalculateERIs_4c_small_molecule(const Eigen::MatrixXd &DMAT) {

          _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
          
          const Eigen::VectorXd& fourc_vector = _fourcenter.get_4c_vector();

          int dftBasisSize = DMAT.rows();
          int vectorSize = (dftBasisSize*(dftBasisSize+1))/2;
          #pragma omp parallel for
          for (int i = 0; i < dftBasisSize; i++) {
            int sum_i = (i*(i+1))/2;
            for (int j = i; j < dftBasisSize; j++) {
              int index_ij = dftBasisSize * i - sum_i + j;
              int index_ij_kl_a = vectorSize * index_ij - (index_ij*(index_ij+1))/2;
              for (int k = 0; k < dftBasisSize; k++) {
                int sum_k = (k*(k+1))/2;
                for (int l = k; l < dftBasisSize; l++) {
                  int index_kl = dftBasisSize * k - sum_k + l;

                  int index_ij_kl = index_ij_kl_a + index_kl;
                  if (index_ij > index_kl) index_ij_kl = vectorSize * index_kl - (index_kl*(index_kl+1))/2 + index_ij;

                  if (l == k) {
                    _ERIs(i, j) += DMAT(k, l) * fourc_vector(index_ij_kl);
                  } else {
                    _ERIs(i, j) += 2. * DMAT(k, l) * fourc_vector(index_ij_kl);
                  }

                }
              }
              _ERIs(j, i) = _ERIs(i, j);
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

          const Eigen::VectorXd& fourc_vector = _fourcenter.get_4c_vector();

          int dftBasisSize = DMAT.rows();
          int vectorSize = (dftBasisSize*(dftBasisSize+1))/2;
          #pragma omp parallel for
          for (int thread = 0; thread < nthreads; ++thread) {
            for (int i = thread; i < dftBasisSize; i+= nthreads) {
              int sum_i = (i*(i+1))/2;
              for (int j = i; j < dftBasisSize; j++) {
                int index_ij = DMAT.cols() * i - sum_i + j;
                int index_ij_kl_a = vectorSize * index_ij - (index_ij*(index_ij+1))/2;
                for (int k = 0; k < dftBasisSize; k++) {
                  int sum_k = (k*(k+1))/2;
                  for (int l = k; l < dftBasisSize; l++) {
                    int index_kl = DMAT.cols() * k - sum_k + l;

                    int _index_ij_kl = index_ij_kl_a + index_kl;
                    if (index_ij > index_kl) _index_ij_kl = vectorSize * index_kl - (index_kl*(index_kl+1))/2 + index_ij;
                    double factorij=1;
                    if(i==j){factorij=0.5;}
                    double factorkl=1;
                    if(l==k){factorkl=0.5;}
                    double factor=factorij*factorkl;
                    EXX_thread[thread](i, l) +=factor*DMAT(j, k) * fourc_vector(_index_ij_kl);
                    EXX_thread[thread](j, l) +=factor*DMAT(i, k) * fourc_vector(_index_ij_kl);
                    EXX_thread[thread](i, k) +=factor*DMAT(j, l) * fourc_vector(_index_ij_kl);
                    EXX_thread[thread](j, k) +=factor*DMAT(i, l) * fourc_vector(_index_ij_kl);
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
        
        
        void ERIs::CalculateERIs_4c_direct(const AOBasis& dftbasis, const Eigen::MatrixXd &DMAT) {

          tensor4d::extent_gen extents;
          
          // Number of shells
          int numShells = dftbasis.getNumofShells();
          
          // Initialize ERIs matrix
          _ERIs = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

          #pragma omp parallel
          { // Begin omp parallel

            Eigen::MatrixXd ERIs_thread = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
            
            #pragma omp for
            for (int iShell_3 = 0; iShell_3 < numShells; iShell_3++) {
              const AOShell& shell_3 = dftbasis.getShell(iShell_3);
              int numFunc_3 = shell_3.getNumFunc();
              for (int iShell_4 = iShell_3; iShell_4 < numShells; iShell_4++) {
                const AOShell& shell_4 = dftbasis.getShell(iShell_4);
                int numFunc_4 = shell_4.getNumFunc();    
                for (int iShell_1 = iShell_3; iShell_1 < numShells; iShell_1++) {
                  const AOShell& shell_1 = dftbasis.getShell(iShell_1);
                  int numFunc_1 = shell_1.getNumFunc();
                  for (int iShell_2 = iShell_1; iShell_2 < numShells; iShell_2++) {
                    const AOShell& shell_2 = dftbasis.getShell(iShell_2);
                    int numFunc_2 = shell_2.getNumFunc();

                    // Pre-screening
                    if (_with_screening && CheckScreen(_screening_eps, shell_1, shell_2, shell_3, shell_4))
                      continue;

                    // Get the current 4c block
                    tensor4d block(extents[ range(0, numFunc_1) ][ range(0, numFunc_2) ][ range(0, numFunc_3)][ range(0, numFunc_4)]);
                    for (int i = 0; i < numFunc_1; ++i) {
                      for (int j = 0; j < numFunc_2; ++j) {
                        for (int k = 0; k < numFunc_3; ++k) {
                          for (int l = 0; l < numFunc_4; ++l) {
                            block[i][j][k][l] = 0.0;
                          }
                        }
                      }
                    }
                    bool nonzero=_fourcenter.FillFourCenterRepBlock(block, shell_1, shell_2, shell_3, shell_4);
                    
                    // If there are only zeros, we don't need to put anything in the ERIs matrix
                    if (!nonzero)
                      continue;

                    // Begin fill ERIs matrix

                    FillERIsBlock(ERIs_thread, DMAT, block, shell_1, shell_2, shell_3, shell_4);

                    // Symmetry 1 <--> 2
                    if (iShell_1 != iShell_2)
                      FillERIsBlock(ERIs_thread, DMAT, block, shell_2, shell_1, shell_3, shell_4);

                    // Symmetry 3 <--> 4
                    if (iShell_3 != iShell_4)
                      FillERIsBlock(ERIs_thread, DMAT, block, shell_1, shell_2, shell_4, shell_3);

                    // Symmetry 1 <--> 2 and 3 <--> 4
                    if (iShell_1 != iShell_2 && iShell_3 != iShell_4)
                      FillERIsBlock(ERIs_thread, DMAT, block, shell_2, shell_1, shell_4, shell_3);

                    // Symmetry (1, 2) <--> (3, 4)
                    if (iShell_1 != iShell_3) {

                      // We need the 'transpose' of block
                      tensor4d block2(extents[ range(0, numFunc_3) ][ range(0, numFunc_4) ][ range(0, numFunc_1)][ range(0, numFunc_2)]);
                      for (int i = 0; i < numFunc_1; ++i) {
                        for (int j = 0; j < numFunc_2; ++j) {
                          for (int k = 0; k < numFunc_3; ++k) {
                            for (int l = 0; l < numFunc_4; ++l) {
                              block2[k][l][i][j] = block[i][j][k][l];
                            }
                          }
                        }
                      }

                      FillERIsBlock(ERIs_thread, DMAT, block2, shell_3, shell_4, shell_1, shell_2);

                      // Symmetry 1 <--> 2
                      if (iShell_1 != iShell_2)
                        FillERIsBlock(ERIs_thread, DMAT, block2, shell_3, shell_4, shell_2, shell_1);

                      // Symmetry 3 <--> 4
                      if (iShell_3 != iShell_4)
                        FillERIsBlock(ERIs_thread, DMAT, block2, shell_4, shell_3, shell_1, shell_2);

                      // Symmetry 1 <--> 2 and 3 <--> 4
                      if (iShell_1 != iShell_2 && iShell_3 != iShell_4)
                        FillERIsBlock(ERIs_thread, DMAT, block2, shell_4, shell_3, shell_2, shell_1);
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
          for (int i = 0; i < DMAT.cols(); i++){
            for (int j = i + 1; j < DMAT.rows(); j++){
              _ERIs(j, i) = _ERIs(i, j);
            }
          }

          CalculateEnergy(DMAT);
          return;
        }
        
        
        void ERIs::FillERIsBlock(Eigen::MatrixXd& ERIsCur, const Eigen::MatrixXd& DMAT,
              const tensor4d& block,
              const AOShell& shell_1, const AOShell& shell_2,
              const AOShell& shell_3, const AOShell& shell_4) {

          for (int iFunc_3 = 0; iFunc_3 < shell_3.getNumFunc(); iFunc_3++) {
            int ind_3 = shell_3.getStartIndex() + iFunc_3;
            for (int iFunc_4 = 0; iFunc_4 < shell_4.getNumFunc(); iFunc_4++) {
              int ind_4 = shell_4.getStartIndex() + iFunc_4;

              // Symmetry
              if (ind_3 > ind_4) {
                continue;
              }

              for (int iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
                int ind_1 = shell_1.getStartIndex() + iFunc_1;
                for (int iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
                  int ind_2 = shell_2.getStartIndex() + iFunc_2;

                  // Symmetry
                  if (ind_1 > ind_2) {
                    continue;
                  }

                  // Symmetry for diagonal elements
                  double multiplier = (ind_1 == ind_2 ? 1.0 : 2.0);
                  // Fill ERIs matrix
                  ERIsCur(ind_3, ind_4) += multiplier * DMAT(ind_1, ind_2) * block[iFunc_1][iFunc_2][iFunc_3][iFunc_4];
                  
                } // End loop over functions in shell 2
              } // End loop over functions in shell 1
            } // End loop over functions in shell 4
          } // End loop over functions in shell 3
          
          return;
        }
        
        
        void ERIs::CalculateERIsDiagonals(const AOBasis& dftbasis) {
          
          tensor4d::extent_gen extents;
          
          // Number of shells
          int numShells = dftbasis.getNumofShells();
          // Total number of functions
          int dftBasisSize = dftbasis.AOBasisSize();
          
          _diagonals = Eigen::MatrixXd::Zero(dftBasisSize,dftBasisSize);
          
          for (int iShell_1 = 0; iShell_1 < numShells; iShell_1++) {
            const AOShell& shell_1 = dftbasis.getShell(iShell_1);
            int numFunc_1 = shell_1.getNumFunc();
            for (int iShell_2 = iShell_1; iShell_2 < numShells; iShell_2++) {
              const AOShell& shell_2 = dftbasis.getShell(iShell_2);
              int numFunc_2 = shell_2.getNumFunc();
              
              // Get the current 4c block
              tensor4d block(extents[ range(0, numFunc_1) ][ range(0, numFunc_2) ][ range(0, numFunc_1)][ range(0, numFunc_2)]);
              for (int i = 0; i < numFunc_1; ++i) {
                for (int j = 0; j < numFunc_2; ++j) {
                  for (int k = 0; k < numFunc_1; ++k) {
                    for (int l = 0; l < numFunc_2; ++l) {
                      block[i][j][k][l] = 0.0;
                    }
                  }
                }
              }
              bool nonzero=_fourcenter.FillFourCenterRepBlock(block, shell_1, shell_2, shell_1, shell_2);
              
              if (!nonzero) {
                continue;
              }

              for (int iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
                int ind_1 = shell_1.getStartIndex() + iFunc_1;
                for (int iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
                  int ind_2 = shell_2.getStartIndex() + iFunc_2;

                  // Symmetry
                  if (ind_1 > ind_2) {
                    continue;
                  }
    
                  _diagonals(ind_1, ind_2) = block[iFunc_1][iFunc_2][iFunc_1][iFunc_2];
                  
                  // Symmetry
                  if (ind_1 != ind_2) {
                    _diagonals(ind_2, ind_1) = _diagonals(ind_1, ind_2);
                  }

                }
              }
            }
          }
          
          return;
        }
        
        
        bool ERIs::CheckScreen(double eps,
            const AOShell& shell_1, const AOShell& shell_2,
            const AOShell& shell_3, const AOShell& shell_4) {
          
          const double eps2=eps * eps;
          
          for (int iFunc_3 = 0; iFunc_3 < shell_3.getNumFunc(); iFunc_3++) {
            int ind_3 = shell_3.getStartIndex() + iFunc_3;
            for (int iFunc_4 = 0; iFunc_4 < shell_4.getNumFunc(); iFunc_4++) {
              int ind_4 = shell_4.getStartIndex() + iFunc_4;

              // Symmetry
              if (ind_3 > ind_4) {
                continue;
              }
              
              for (int iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
                int ind_1 = shell_1.getStartIndex() + iFunc_1;
                for (int iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
                  int ind_2 = shell_2.getStartIndex() + iFunc_2;

                  // Symmetry
                  if (ind_1 > ind_2) {
                    continue;
                  }
                  
                  // Cauchy–Schwarz
                  // <ab|cd> <= sqrt(<ab|ab>) * sqrt(<cd|cd>)
                  double ub = _diagonals(ind_1, ind_2) * _diagonals(ind_3, ind_4);
                  
                  // Compare with tolerance
                  if (ub > eps2) {
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
          _EXXsenergy=_EXXs.cwiseProduct(DMAT).sum();
          return;
        }
        
        
        void ERIs::printERIs(){
          for (int i=0; i< _ERIs.cols(); i++){
            for (int j=0; j< _ERIs.rows();j++){
              std::cout << "ERIs [" << i<<":"<<j<<"]="<<_ERIs(i,j)<<std::endl;
            }
          }
          return;
        }
        
    }
}