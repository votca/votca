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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/bse.h>
#include <votca/tools/linalg.h>
using boost::format;
using std::flush;

namespace votca {
  namespace xtp {

    void BSE::Solve_triplets() {

      // add full QP Hamiltonian contributions to free transitions
      MatrixXfd H = MatrixXfd::Zero(_bse_size,_bse_size);
      Add_Hd<real_gwbse>(H);
      Add_Hqp<real_gwbse>(H);
      
      CTP_LOG(ctp::logDEBUG, *_log)
        << ctp::TimeStamp() << " Setup TDA triplet hamiltonian " << flush;
      CTP_LOG(ctp::logDEBUG, *_log)
        << ctp::TimeStamp() << " Solving for first "<<_bse_nmax<<" eigenvectors"<< flush;
      tools::linalg_eigenvalues(H , _bse_triplet_energies, _bse_triplet_coefficients ,_bse_nmax );
      return;
    }

    void BSE::Solve_singlets() {
       
      MatrixXfd H = MatrixXfd::Zero(_bse_size,_bse_size);
      Add_Hd<real_gwbse>(H);
      Add_Hqp<real_gwbse>(H);
      Add_Hx<real_gwbse>(H,2.0);
      CTP_LOG(ctp::logDEBUG, *_log)
        << ctp::TimeStamp() << " Setup TDA singlet hamiltonian " << flush;
      CTP_LOG(ctp::logDEBUG, *_log)
        << ctp::TimeStamp() << " Solving for first "<<_bse_nmax<<" eigenvectors"<< flush;
      tools::linalg_eigenvalues(H, _bse_singlet_energies, _bse_singlet_coefficients , _bse_nmax );
      return;
    }
    
     void BSE::SetupHs(){
      _eh_s = MatrixXfd::Zero(_bse_size,_bse_size);
      Add_Hd<real_gwbse>(_eh_s);
      Add_Hqp<real_gwbse>(_eh_s);
      Add_Hx<real_gwbse>(_eh_s,2.0);
     }
  
  void BSE::SetupHt(){
      _eh_t = MatrixXfd::Zero(_bse_size,_bse_size);
      Add_Hd<real_gwbse>(_eh_t);
      Add_Hqp<real_gwbse>(_eh_t);
  }
    

    void BSE::Solve_singlets_BTDA() {

      // For details of the method, see EPL,78(2007)12001,
      // Nuclear Physics A146(1970)449, Nuclear Physics A163(1971)257.

      // setup resonant (A) and RARC blocks (B)
        
       //corresponds to 
      // _ApB = (_eh_d +_eh_qp + _eh_d2 + 4.0 * _eh_x);
      // _AmB = (_eh_d +_eh_qp - _eh_d2);
        Eigen::MatrixXd _ApB=Eigen::MatrixXd::Zero(_bse_size,_bse_size);
        Add_Hd<double>(_ApB);
        Add_Hqp<double>(_ApB);
        
        
        Eigen::MatrixXd _AmB=_ApB;
        Add_Hd2<double>(_AmB,-1.0);
        
        Add_Hx<double>(_ApB,4.0);
        Add_Hd2<double>(_ApB,1.0);
        CTP_LOG(ctp::logDEBUG, *_log)
        << ctp::TimeStamp() << " Setup singlet hamiltonian " << flush;
     

      // calculate Cholesky decomposition of A-B = LL^T. It throws an error if not positive definite
      //(A-B) is not needed any longer and can be overwritten
      CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Trying Cholesky decomposition of KAA-KAB" << flush;
      Eigen::LLT< Eigen::Ref<Eigen::MatrixXd> > L(_AmB);
      
       for (int i=0;i<_AmB.rows();++i){
          for (int j=i+1;j<_AmB.cols();++j){
          _AmB(i,j)=0;
          }
        }
            
      std::string success="successful";
      if(L.info()!=0){
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() <<" Cholesky decomposition of KAA-KAB was unsucessful. Try a smaller basisset. This can indicate a triplet instability."<<flush;
        throw std::runtime_error("Cholesky decompostion failed");
      }else{
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() <<" Cholesky decomposition of KAA-KAB was successful"<<flush;
      }
      
      Eigen::MatrixXd temp= _ApB*_AmB;
      _ApB.noalias() =_AmB.transpose() *temp;
      temp.resize(0,0);
      CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculated H = L^T(A+B)L " << flush;
      Eigen::VectorXd eigenvalues;
      Eigen::MatrixXd eigenvectors;
      
      CTP_LOG(ctp::logDEBUG, *_log)
        << ctp::TimeStamp() << " Solving for first "<<_bse_nmax<<" eigenvectors"<< flush;
      bool success_diag=tools::linalg_eigenvalues(_ApB, eigenvalues, eigenvectors ,_bse_nmax);
      if(!success_diag){
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Could not solve problem" << flush;
      }else{
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;
      }
      _ApB.resize(0,0);
      eigenvalues=eigenvalues.cwiseSqrt();
     
      #if (GWBSE_DOUBLE)
      _bse_singlet_energies =eigenvalues;
#else
      _bse_singlet_energies = eigenvalues.cast<float>(); 
#endif

      // reconstruct real eigenvectors X_l = 1/2 [sqrt(eps_l) (L^T)^-1 + 1/sqrt(eps_l)L ] R_l
      //                               Y_l = 1/2 [sqrt(eps_l) (L^T)^-1 - 1/sqrt(eps_l)L ] R_l

      // determine inverse of L^T
       
     Eigen::MatrixXd LmT = _AmB.inverse().transpose();

      int dim = LmT.rows();
      _bse_singlet_energies.resize(_bse_nmax);
      _bse_singlet_coefficients.resize(dim, _bse_nmax); // resonant part (_X_evec)
      _bse_singlet_coefficients_AR.resize(dim, _bse_nmax); // anti-resonant part (_Y_evec)
      for (int _i = 0; _i < _bse_nmax; _i++) {

        //real_gwbse sqrt_eval = sqrt(_eigenvalues(_i));
        double sqrt_eval = sqrt(_bse_singlet_energies(_i));
        // get l-th reduced EV
          #if (GWBSE_DOUBLE)
        _bse_singlet_coefficients.col(_i) = (0.5 / sqrt_eval * (_bse_singlet_energies(_i) * LmT + _AmB) * eigenvectors.col(_i));
        _bse_singlet_coefficients_AR.col(_i) = (0.5 / sqrt_eval * (_bse_singlet_energies(_i) * LmT - _AmB) * eigenvectors.col(_i));
#else
        _bse_singlet_coefficients.col(_i) = (0.5 / sqrt_eval * (_bse_singlet_energies(_i) * LmT + _AmB) * eigenvectors.col(_i)).cast<float>();
        _bse_singlet_coefficients_AR.col(_i) = (0.5 / sqrt_eval * (_bse_singlet_energies(_i) * LmT - _AmB) * eigenvectors.col(_i)).cast<float>();
#endif

      }

      return;
    }

   
template <typename T>
    void BSE::Add_Hqp(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {
    
    const Eigen::MatrixXd& Hqp=*_Hqp;

#pragma omp parallel for
      for (size_t _v1 = 0; _v1 < _bse_vtotal; _v1++) {
        for (size_t _c1 = 0; _c1 < _bse_ctotal; _c1++) {
          size_t _index_vc = _bse_ctotal * _v1 + _c1;
          // diagonal
          H(_index_vc, _index_vc) += Hqp(_c1 + _bse_vtotal, _c1 + _bse_vtotal) -Hqp(_v1, _v1);
          // v->c
          for (size_t _c2 = 0; _c2 < _bse_ctotal; _c2++) {
            size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
            if (_c1 != _c2) {
              H(_index_vc, _index_vc2) += Hqp(_c1 + _bse_vtotal, _c2 + _bse_vtotal);
            }
          }

          // c-> v
          for (size_t _v2 = 0; _v2 < _bse_vtotal; _v2++) {
            size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
            if (_v1 != _v2) {
              H(_index_vc, _index_vc2) -= Hqp(_v1, _v2);
            }
          }
        }
      }
      return;
    }

template <typename T>
    void BSE::Add_Hd(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {
      // gwbasis size
      size_t auxsize = _Mmn->getAuxDimension();
      size_t vxv_size=_bse_vtotal * _bse_vtotal;
      size_t cxc_size=_bse_ctotal * _bse_ctotal;

      
      // messy procedure, first get two matrices for occ and empty subbparts
      // store occs directly transposed
      MatrixXfd _storage_v = MatrixXfd::Zero(auxsize, vxv_size);
#pragma omp parallel for
      for (size_t _v1 = 0; _v1 < _bse_vtotal; _v1++) {
        const MatrixXfd& Mmn = (*_Mmn)[_v1 + _bse_vmin ];
        for (size_t _i_gw = 0; _i_gw < auxsize; _i_gw++) {
          for (size_t _v2 = 0; _v2 < _bse_vtotal; _v2++) {
            size_t _index_vv = _bse_vtotal * _v1 + _v2;         
              _storage_v(_i_gw,_index_vv) = Mmn( _v2 + _bse_vmin,_i_gw);
            }
        }
      }

      MatrixXfd _storage_c = MatrixXfd::Zero(auxsize,cxc_size);
#pragma omp parallel for
      for (size_t _c1 = 0; _c1 < _bse_ctotal; _c1++) {
        const MatrixXfd& Mmn = (*_Mmn)[_c1 + _bse_cmin];
        for (size_t _i_gw = 0; _i_gw < auxsize; _i_gw++) {
          for (size_t _c2 = 0; _c2 < _bse_ctotal; _c2++) {
            size_t _index_cc = _bse_ctotal * _c1 + _c2;
            _storage_c(_i_gw, _index_cc) = Mmn( _c2 + _bse_cmin,_i_gw);
          }
        }
      }

      // store elements in a vtotal^2 x ctotal^2 matrix
      MatrixXfd _storage_prod = _storage_v.transpose() *_storage_c;

      // now patch up _storage for screened interaction
#pragma omp parallel for
      for (size_t _i_gw = 0; _i_gw < auxsize; _i_gw++) {
        if (_ppm->getPpm_weight()(_i_gw) < 1.e-9) {
          for (size_t _v = 0; _v < vxv_size; _v++) {
            _storage_v(_i_gw,_v ) = 0;
          }
          for (size_t _c = 0; _c < cxc_size; _c++) {
            _storage_c(_i_gw, _c) = 0;
          }

        } else {
          double ppm_factor = sqrt(_ppm->getPpm_weight()(_i_gw));
          for (size_t _v = 0; _v < vxv_size; _v++) {
            _storage_v(_i_gw,_v ) *= ppm_factor;
          }
          for (size_t _c = 0; _c < cxc_size; _c++) {
            _storage_c(_i_gw, _c) *= ppm_factor;
          }
        }
      }

      // multiply and subtract from _storage_prod

      _storage_prod -= _storage_v.transpose()*_storage_c;

      // free storage_v and storage_c
      _storage_c.resize(0, 0);
      _storage_v.resize(0, 0);

      // finally resort into _eh_d
      
#pragma omp parallel for
      for (size_t _v1 = 0; _v1 < _bse_vtotal; _v1++) {
        for (size_t _v2 = 0; _v2 < _bse_vtotal; _v2++) {
          size_t _index_vv = _bse_vtotal * _v1 + _v2;

          for (size_t _c1 = 0; _c1 < _bse_ctotal; _c1++) {
            size_t _index_vc1 = _bse_ctotal * _v1 + _c1;

            for (size_t _c2 = 0; _c2 < _bse_ctotal; _c2++) {
              size_t _index_vc2 = _bse_ctotal * _v2 + _c2;
              size_t _index_cc = _bse_ctotal * _c1 + _c2;

              H(_index_vc1, _index_vc2) -= _storage_prod(_index_vv, _index_cc);
            }
          }
        }
      }

      return;
    }
template <typename T>
    void BSE::Add_Hd2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H, double factor) {
      // gwbasis size
      size_t auxsize = _Mmn->getAuxDimension();
      size_t bse_vxc_total=_bse_vtotal * _bse_ctotal;
      // messy procedure, first get two matrices for occ and empty subbparts
      // store occs directly transposed
      MatrixXfd _storage_cv = MatrixXfd::Zero(auxsize,bse_vxc_total);
#pragma omp parallel for
      for (size_t _c1 = 0; _c1 < _bse_ctotal; _c1++) {
        const MatrixXfd& Mmn = (*_Mmn)[_c1 + _bse_cmin ];
        for (size_t _i_gw = 0; _i_gw < auxsize; _i_gw++) {
          for (size_t _v2 = 0; _v2 < _bse_vtotal; _v2++) {
            size_t _index_cv = _bse_vtotal * _c1 + _v2;
            _storage_cv(_i_gw,_index_cv ) = Mmn( _v2 + _bse_vmin,_i_gw);
          }
        }
      }

      MatrixXfd _storage_vc = MatrixXfd::Zero(auxsize, bse_vxc_total);
#pragma omp parallel for
      for (size_t _v1 = 0; _v1 < _bse_vtotal; _v1++) {
        const MatrixXfd& Mmn = (*_Mmn)[_v1 + _bse_vmin];
        for (size_t _i_gw = 0; _i_gw < auxsize; _i_gw++) {
          for (size_t _c2 = 0; _c2 < _bse_ctotal; _c2++) {
            size_t _index_vc = _bse_ctotal * _v1 + _c2;
            _storage_vc(_i_gw, _index_vc) = Mmn(_c2 + _bse_cmin,_i_gw);
          }
        }
      }

      // store elements in a vtotal^2 x ctotal^2 matrix
      MatrixXfd _storage_prod = _storage_cv.transpose()* _storage_vc;


      // now patch up _storage for screened interaction
#pragma omp parallel for
      for (size_t _i_gw = 0; _i_gw < auxsize; _i_gw++) {
        if (_ppm->getPpm_weight()(_i_gw) < 1.e-9) {
          for (size_t _v = 0; _v < bse_vxc_total; _v++) {
            _storage_vc(_i_gw, _v) = 0;
          }
          for (size_t _c = 0; _c < bse_vxc_total; _c++) {
            _storage_cv(_i_gw,_c ) = 0;
          }
        } else {
          double ppm_factor = sqrt(_ppm->getPpm_weight()(_i_gw));
          for (size_t _v = 0; _v < bse_vxc_total; _v++) {
            _storage_vc(_i_gw, _v) *= ppm_factor;
          }
          for (size_t _c = 0; _c < bse_vxc_total; _c++) {
            _storage_cv(_i_gw,_c ) *= ppm_factor;
          }
        }
      }

      // multiply and subtract from _storage_prod
      _storage_prod -= _storage_cv.transpose() * _storage_vc;

      // free storage_v and storage_c
      _storage_cv.resize(0, 0);
      _storage_vc.resize(0, 0);
  
#pragma omp parallel for
      for (size_t _v1 = 0; _v1 < _bse_vtotal; _v1++) {
        for (size_t _v2 = 0; _v2 < _bse_vtotal; _v2++) {
          for (size_t _c1 = 0; _c1 < _bse_ctotal; _c1++) {
            size_t _index_v1c1 = _bse_ctotal * _v1 + _c1;

            size_t _index_c1v2 = _bse_vtotal * _c1 + _v2;

            for (size_t _c2 = 0; _c2 < _bse_ctotal; _c2++) {
              size_t _index_v2c2 = _bse_ctotal * _v2 + _c2;
              size_t _index_v1c2 = _bse_ctotal * _v1 + _c2;

              H(_index_v1c1, _index_v2c2) -= factor*_storage_prod(_index_c1v2, _index_v1c2);

            }
          }
        }
      }
      return;
    }
template <typename T>
    void BSE::Add_Hx(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H, double factor) {

    
      // gwbasis size
      size_t auxsize = _Mmn->getAuxDimension();

      // get a different storage for 3-center integrals we need
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _storage = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(auxsize, _bse_size);

      // occupied levels
#pragma omp parallel for
      for (size_t _v = 0; _v < _bse_vtotal; _v++) {

        const MatrixXfd& Mmn = (*_Mmn)[_v + _bse_vmin];
        // empty levels
        for (size_t _i_gw = 0; _i_gw < auxsize; _i_gw++) {
          for (size_t _c = 0; _c < _bse_ctotal; _c++) {
            size_t _index_vc = _bse_ctotal * _v + _c;
            _storage(_i_gw, _index_vc) = Mmn(_c + _bse_cmin,_i_gw);
          }
        }
      }

      // with this storage, _eh_x is obtained by matrix multiplication
      H += factor*_storage.transpose() * _storage;
      return;
    }

    void BSE::printFragInfo(Population& pop, int i){
      CTP_LOG(ctp::logINFO, *_log) << (format("           Fragment A -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
              % (100.0 * pop.popH[i](0)) % (100.0 * pop.popE[i](0)) % (pop.Crgs[i](0)) % (pop.Crgs[i](0) + pop.popGs(0))).str() << flush;
      CTP_LOG(ctp::logINFO, *_log) << (format("           Fragment B -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
              % (100.0 * pop.popH[i](1)) % (100.0 * pop.popE[i](1)) % (pop.Crgs[i](1)) % (pop.Crgs[i](1) + pop.popGs(1))).str() << flush;
      return;
    }

    void BSE::printWeights(unsigned i_bse, double weight){
      if (weight > _min_print_weight) {
        CTP_LOG(ctp::logINFO, *_log) << (format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%")
                % (_homo - _index2v[i_bse]) % (_index2c[i_bse] - _homo - 1) % (100.0 * weight)).str() << flush;
      }
      return;
    }

    void BSE::Analyze_singlets(const AOBasis& dftbasis) {

      Interaction act;
      Population pop;
      
      std::vector< tools::vec > transition_dipoles=CalcCoupledTransition_Dipoles(dftbasis);
      _orbitals.TransitionDipoles()=transition_dipoles;
      std::vector<double> oscs = _orbitals.Oscillatorstrengths();
      
      if(tools::globals::verbose){
        act = Analyze_eh_interaction("singlet");
      }
      if(dftbasis.getAOBasisFragA()>0){
        pop=FragmentPopulations("singlet",dftbasis);
        _orbitals.setFragmentChargesSingEXC(pop.Crgs);
        _orbitals.setFragment_E_localisation_singlet(pop.popE);
        _orbitals.setFragment_H_localisation_singlet(pop.popH);
        _orbitals.setFragmentChargesGS(pop.popGs);
      }
      
      double hrt2ev = tools::conv::hrt2ev;
      CTP_LOG(ctp::logINFO, *_log) << (format("  ====== singlet energies (eV) ====== ")).str() << flush;
      int maxoutput=(_bse_nmax>200) ? 200:_bse_nmax;
      for (int i = 0; i < maxoutput; ++i) {
       
        const tools::vec& trdip = transition_dipoles[i];
        double osc = oscs[i];
        if (tools::globals::verbose) {
          CTP_LOG(ctp::logINFO, *_log) << (format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_x> = %5$+1.4f <K_d> = %6$+1.4f")
                  % (i + 1) % (hrt2ev * _bse_singlet_energies(i)) % (1240.0 / (hrt2ev * _bse_singlet_energies(i)))
                  % (hrt2ev * act.qp_contrib(i)) % (hrt2ev * act.exchange_contrib(i)) % (hrt2ev * act.direct_contrib(i))).str() << flush;
        } else {
          CTP_LOG(ctp::logINFO, *_log) << (format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                  % (i + 1) % (hrt2ev * _bse_singlet_energies(i)) % (1240.0 / (hrt2ev * _bse_singlet_energies(i)))).str() << flush;
        }

        CTP_LOG(ctp::logINFO, *_log) << (format("           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f")
                % trdip.getX() % trdip.getY() % trdip.getZ() % (trdip * trdip) % osc).str() << flush;
        for (unsigned i_bse = 0; i_bse < _bse_size; ++i_bse) {
          // if contribution is larger than 0.2, print
          double weight = pow(_bse_singlet_coefficients(i_bse, i), 2);
          if (_bse_singlet_coefficients_AR.rows()>0){
               weight-= pow(_bse_singlet_coefficients_AR(i_bse, i), 2);
          }
          printWeights(i_bse, weight);
        }
        // results of fragment population analysis 
        if (dftbasis.getAOBasisFragA() > 0 && dftbasis.getAOBasisFragB()>0) {
          printFragInfo(pop, i);
        }

        CTP_LOG(ctp::logINFO, *_log) << flush;
      }
      return;
    }
    
    
    
    
    void BSE::Analyze_triplets(const AOBasis& dftbasis) {

      Interaction act;
      Population pop;
            
      if(tools::globals::verbose){
        act = Analyze_eh_interaction("triplet");
      }
      if(dftbasis.getAOBasisFragA()>0){
        pop=FragmentPopulations("triplet",dftbasis);
        _orbitals.setFragmentChargesTripEXC(pop.Crgs);
        _orbitals.setFragment_E_localisation_triplet(pop.popE);
        _orbitals.setFragment_H_localisation_triplet(pop.popH);
        _orbitals.setFragmentChargesGS(pop.popGs);
      }
      CTP_LOG(ctp::logINFO, *_log) << (format("  ====== triplet energies (eV) ====== ")).str() << flush;
      int maxoutput=(_bse_nmax>200) ? 200:_bse_nmax;
      for (int i = 0; i < maxoutput; ++i) {
        
        if (tools::globals::verbose) {
          CTP_LOG(ctp::logINFO, *_log) << (format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_d> = %5$+1.4f")
                  % (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) % (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i)))
                  % (tools::conv::hrt2ev * act.qp_contrib(i)) % (tools::conv::hrt2ev *act.direct_contrib(i))).str() << flush;
        } else {
          CTP_LOG(ctp::logINFO, *_log) << (format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                  % (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) % (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i)))).str() << flush;
        }

        for (unsigned i_bse = 0; i_bse < _bse_size; ++i_bse) {
          // if contribution is larger than 0.2, print
          double weight = pow(_bse_triplet_coefficients(i_bse, i), 2);
          printWeights(i_bse, weight);
        }
        // results of fragment population analysis 
        if (dftbasis.getAOBasisFragA() > 0) {
          printFragInfo(pop, i);
        }
        CTP_LOG(ctp::logINFO, *_log) << (format("   ")).str() << flush;
      }

      // storage to orbitals object

      return;
    }

    Eigen::VectorXd BSE::Analyze_IndividualContribution(const std::string& spin,const MatrixXfd& H){
        Eigen::VectorXd contrib=Eigen::VectorXd::Zero(_bse_nmax);
        if (spin == "singlet") {
            for (int i_exc = 0; i_exc < _bse_nmax; i_exc++) {
                MatrixXfd _slice_R = _bse_singlet_coefficients.block(0, i_exc, _bse_size, 1);
                contrib(i_exc) =  (_slice_R.transpose()*H * _slice_R).value();
                if (_bse_singlet_coefficients_AR.cols() > 0) {
                    MatrixXfd _slice_AR = _bse_singlet_coefficients_AR.block(0, i_exc, _bse_size, 1);
                    // get anti-resonant contribution from direct Keh 
                    contrib(i_exc)-= (_slice_AR.transpose()*H * _slice_AR).value();           
                }
            }
        } else if (spin == "triplet") {
            for (int i_exc = 0; i_exc < _bse_nmax; i_exc++) {
                MatrixXfd _slice_R = _bse_triplet_coefficients.block(0, i_exc, _bse_size, 1);
                contrib(i_exc) =  (_slice_R.transpose()*H * _slice_R).value();
            }
        } else {
            throw std::runtime_error("BSE::Analyze_eh_interaction:Spin not known!");
        }
        return contrib;
    }

    BSE::Interaction BSE::Analyze_eh_interaction(const std::string& spin) {

      Interaction analysis;
      MatrixXfd H = MatrixXfd::Zero(_bse_size, _bse_size);
      Add_Hqp(H); 
      analysis.qp_contrib=Analyze_IndividualContribution(spin,H);
      
      H = MatrixXfd::Zero(_bse_size, _bse_size);
      Add_Hd(H);
      analysis.direct_contrib=Analyze_IndividualContribution(spin,H);
      if (spin == "singlet") {
          H = MatrixXfd::Zero(_bse_size, _bse_size);
          Add_Hx(H,2.0);
          analysis.exchange_contrib=Analyze_IndividualContribution(spin,H);
      }else{
            analysis.exchange_contrib=Eigen::VectorXd::Zero(0);
      }
      
      return analysis;
    }

    BSE::Population BSE::FragmentPopulations(const std::string& spin, const AOBasis& dftbasis) {
      Population pop;
      // Mulliken fragment population analysis
        // get overlap matrix for DFT basisset
        AOOverlap _dftoverlap;
        // Fill overlap
        _dftoverlap.Fill(dftbasis);
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftoverlap.Matrix().rows() << flush;
        // ground state populations
        Eigen::MatrixXd DMAT = _orbitals.DensityMatrixGroundState();

        Eigen::VectorXd nuccharges = _orbitals.FragmentNuclearCharges(dftbasis.getAOBasisFragA());
        Eigen::VectorXd pops = _orbitals.LoewdinPopulation(DMAT, _dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
        pop.popGs=nuccharges - pops;
        // population to electron charges and add nuclear charges         
        for (int _i_state = 0; _i_state < _bse_nmax; _i_state++) {

          // checking Density Matrices
          std::vector< Eigen::MatrixXd > DMAT = _orbitals.DensityMatrixExcitedState(spin, _i_state);
          // hole part
          Eigen::VectorXd popsH = _orbitals.LoewdinPopulation(DMAT[0], _dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
          pop.popH.push_back(popsH);
          // electron part
          Eigen::VectorXd popsE = _orbitals.LoewdinPopulation(DMAT[1], _dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
          pop.popE.push_back(popsE);
          // update effective charges
          Eigen::VectorXd diff = popsH - popsE;
          pop.Crgs.push_back(diff);
        }
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Ran Excitation fragment population analysis " << flush;
     
      return pop;
    }

    std::vector<Eigen::MatrixXd > BSE::CalcFreeTransition_Dipoles(const AOBasis& dftbasis) {
      const Eigen::MatrixXd& dft_orbitals = _orbitals.MOCoefficients();
      // Testing electric dipole AOMatrix
      AODipole _dft_dipole;
      _dft_dipole.Fill(dftbasis);

      // now transition dipole elements for free interlevel transitions
      std::vector<Eigen::MatrixXd > interlevel_dipoles;

      Eigen::MatrixXd empty = dft_orbitals.block(0,_bse_cmin,dftbasis.AOBasisSize() , _bse_ctotal);
      Eigen::MatrixXd occ = dft_orbitals.block(0,_bse_vmin, dftbasis.AOBasisSize(), _bse_vtotal);
      for (int _i_comp = 0; _i_comp < 3; _i_comp++) {
        interlevel_dipoles.push_back(occ.transpose() * _dft_dipole.Matrix()[_i_comp] * empty);
      }
      CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculated free interlevel transition dipole moments " << flush;
      return interlevel_dipoles;
    }

    std::vector<tools::vec > BSE::CalcCoupledTransition_Dipoles(const AOBasis& dftbasis) {
    std::vector<Eigen::MatrixXd > interlevel_dipoles= CalcFreeTransition_Dipoles(dftbasis);
    std::vector<tools::vec > dipols;
      double sqrt2 = sqrt(2.0);
      for (int _i_exc = 0; _i_exc < _bse_nmax; _i_exc++) {
        tools::vec _tdipole = tools::vec(0, 0, 0);

        for (unsigned _v = 0; _v < _bse_vtotal; _v++) {
          for (unsigned _c = 0; _c < _bse_ctotal; _c++) {
            int index_vc = _bse_ctotal * _v + _c;
            double factor = sqrt2 * _bse_singlet_coefficients(index_vc, _i_exc);
            if (_bse_singlet_coefficients_AR.rows()>0) {
              factor += sqrt2 * _bse_singlet_coefficients_AR(index_vc, _i_exc);
            }
            // The Transition dipole is sqrt2 bigger because of the spin, the excited state is a linear combination of 2 slater determinants, where either alpha or beta spin electron is excited
            _tdipole.x() += factor * interlevel_dipoles[0](_v, _c);
            _tdipole.y() += factor * interlevel_dipoles[1](_v, _c);
            _tdipole.z() += factor * interlevel_dipoles[2](_v, _c);
          }

        }
        dipols.push_back(_tdipole);
       
      }
      return dipols;
    }
    


  }
};
