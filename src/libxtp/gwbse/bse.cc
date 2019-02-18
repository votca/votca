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

#include <votca/xtp/davidsonsolver.h>
#include <votca/xtp/bse_operator.h>

#include "votca/xtp/qmstate.h"
#include "votca/xtp/vc2index.h"

using boost::format;
using std::flush;

namespace votca {
  namespace xtp {
  
    void BSE::Solve_triplets() {

      BSE_OPERATOR Ht(_orbitals, _log, _Mmn, _Hqp);
      BSE::configure_operator(Ht);

      Ht.setHqp(1);
      Ht.setHd(1);
      
      CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Setup TDA triplet hamiltonian " << flush;

      Eigen::VectorXd _tmp_es;
      Eigen::MatrixXd _tmp_mo;
      BSE::solve_hermitian(Ht, _tmp_es, _tmp_mo);

      #if (GWBSE_DOUBLE)
        _bse_triplet_energies = _tmp_es;
        _bse_triplet_coefficients = _tmp_mo;  
      #else
        _bse_triplet_energies = _tmp_es.cast<float>();
        _bse_triplet_coefficients = _tmp_mo.cast<float>();  
      #endif

      return;
    }

    void BSE::Solve_singlets() {
        if(_opt.useTDA){
            Solve_singlets_TDA();
        }
        else {
          if (_opt.davidson){
            CTP_LOG(ctp::logDEBUG, _log)
              << ctp::TimeStamp() << " Davidson solver not implemented for BTDA. Full diagonalization." << flush;
            _opt.davidson = 0;
          }
          Solve_singlets_BTDA(); 
        }
    }

    void BSE::Solve_singlets_TDA() {

      BSE_OPERATOR Hs(_orbitals, _log, _Mmn, _Hqp);
      BSE::configure_operator(Hs);

      Hs.setHx(2);
      Hs.setHqp(1);
      Hs.setHd(1);

      CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Setup TDA singlet hamiltonian " << flush;

      Eigen::VectorXd _tmp_es;
      Eigen::MatrixXd _tmp_mo;
      BSE::solve_hermitian(Hs, _tmp_es, _tmp_mo);

      #if (GWBSE_DOUBLE)
        _bse_singlet_energies = _tmp_es;
        _bse_singlet_coefficients = _tmp_mo;  
      #else
        _bse_singlet_energies = _tmp_es.cast<float>();
        _bse_singlet_coefficients = _tmp_mo.cast<float>();  
      #endif
      
    }

    void BSE::SetupHs() {
      BSE_OPERATOR Hs(_orbitals, _log, _Mmn, _Hqp);
      BSE::configure_operator(Hs);
      

      Hs.setHx(2);
      Hs.setHqp(1);
      Hs.setHd(1); 
      _eh_s = Hs.get_full_matrix();
    }
    
    void BSE::SetupHt() {
      BSE_OPERATOR Ht(_orbitals, _log, _Mmn, _Hqp);
      BSE::configure_operator(Ht);
      

      Ht.setHqp(1);
      Ht.setHd(1); 
      _eh_t = Ht.get_full_matrix();
    }

    void BSE::configure_operator(BSE_OPERATOR &h) {
      h._opt.homo = this->_opt.homo;
      h._opt.rpamin = this->_opt.rpamin;
      h._opt.rpamax = this->_opt.rpamax;
      h._opt.qpmin = this->_opt.qpmin;
      h._opt.vmin = this->_opt.vmin;
      h._opt.cmax = this->_opt.cmax;
      h._bse_cmin = _opt.homo+1;
      h._bse_vtotal = _bse_vmax - _opt.vmin + 1;
      h._bse_ctotal =_opt.cmax - _bse_cmin + 1;
      h._bse_size = _bse_vtotal * _bse_ctotal; 
      h._size = _bse_size;
      h.SetupDirectInteractionOperator();
      return ;
    }

    void BSE::solve_hermitian(BSE_OPERATOR &h, Eigen::VectorXd &energies, Eigen::MatrixXd &coefficients) {

      CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Solving for first "<<_opt.nmax<<" eigenvectors"<< flush;

      if (_opt.davidson) {
        
        DavidsonSolver DS(_log);
        DS.set_correction(_opt.davidson_correction);
      
        if(_opt.matrixfree) {
          CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Using matrix free method"<< flush;
          DS.solve(h,_opt.nmax);
        } else {
          CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Using full matrix method"<< flush;
          Eigen::MatrixXd hfull = h.get_full_matrix();
          DS.solve(hfull,_opt.nmax);
        }
        energies = DS.eigenvalues(); 
        coefficients = DS.eigenvectors();
      }

      else {
        CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Lapack Diagonalization" << flush;
        Eigen::MatrixXd hfull = h.get_full_matrix();
        tools::linalg_eigenvalues(hfull, energies, coefficients , _opt.nmax );
      }
      return;
    }

    void BSE::Solve_singlets_BTDA() {

      // For details of the method, see EPL,78(2007)12001,
      // Nuclear Physics A146(1970)449, Nuclear Physics A163(1971)257.
      // setup resonant (A) and RARC blocks (B)
       //corresponds to 
      // _ApB = (_eh_d +_eh_qp + _eh_d2 + 4.0 * _eh_x);
      // _AmB = (_eh_d +_eh_qp - _eh_d2);

      BSE_OPERATOR Hs_ApB(_orbitals, _log, _Mmn, _Hqp);
      BSE::configure_operator(Hs_ApB);
      

      Hs_ApB.setHd(1);
      Hs_ApB.setHqp(1);
      Hs_ApB.setHd2(1);
      Hs_ApB.setHx(4);

      BSE_OPERATOR Hs_AmB(_orbitals, _log, _Mmn, _Hqp);
      BSE::configure_operator(Hs_AmB);
      
      Hs_AmB.setHd(1);
      Hs_AmB.setHqp(1);
      Hs_AmB.setHd2(-1);

      Eigen::MatrixXd ApB = Hs_ApB.get_full_matrix();
      Eigen::MatrixXd AmB = Hs_AmB.get_full_matrix();

      CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Setup singlet hamiltonian " << flush;
     
      // calculate Cholesky decomposition of A-B = LL^T. It throws an error if not positive definite
      //(A-B) is not needed any longer and can be overwritten
      CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Trying Cholesky decomposition of KAA-KAB" << flush;
      Eigen::LLT< Eigen::Ref<Eigen::MatrixXd> > L(AmB);
      
       for (int i=0;i<AmB.rows();++i){
          for (int j=i+1;j<AmB.cols();++j){
          AmB(i,j)=0;
          }
        }

      if(L.info()!=0){
        CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() <<" Cholesky decomposition of KAA-KAB was unsucessful. Try a smaller basisset. This can indicate a triplet instability."<<flush;
        throw std::runtime_error("Cholesky decompostion failed");
      }else{
        CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() <<" Cholesky decomposition of KAA-KAB was successful"<<flush;
      }
      
      Eigen::MatrixXd temp= ApB*AmB;
      ApB.noalias() =AmB.transpose()*temp;
      temp.resize(0,0);
      CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Calculated H = L^T(A+B)L " << flush;
      Eigen::VectorXd eigenvalues;
      Eigen::MatrixXd eigenvectors;     
      CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Solving for first "<<_opt.nmax<<" eigenvectors"<< flush;
      bool success_diag=tools::linalg_eigenvalues(ApB, eigenvalues, eigenvectors ,_opt.nmax);
      if(!success_diag){
        CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Could not solve problem" << flush;
      }else{
        CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;
      }
      ApB.resize(0,0);
      eigenvalues=eigenvalues.cwiseSqrt();
     
#if (GWBSE_DOUBLE)
      _bse_singlet_energies = eigenvalues;
#else
      _bse_singlet_energies = eigenvalues.cast<float>(); 
#endif
      // reconstruct real eigenvectors X_l = 1/2 [sqrt(eps_l) (L^T)^-1 + 1/sqrt(eps_l)L ] R_l
      //                               Y_l = 1/2 [sqrt(eps_l) (L^T)^-1 - 1/sqrt(eps_l)L ] R_l
      // determine inverse of L^T
     Eigen::MatrixXd LmT = AmB.inverse().transpose();
      int dim = LmT.rows();
      _bse_singlet_energies.resize(_opt.nmax);
      _bse_singlet_coefficients.resize(dim, _opt.nmax); // resonant part (_X_evec)
      _bse_singlet_coefficients_AR.resize(dim, _opt.nmax); // anti-resonant part (_Y_evec)
      for (int level = 0; level < _opt.nmax; level++) {
        double sqrt_eval = std::sqrt(_bse_singlet_energies(level));
        // get l-th reduced EV
#if (GWBSE_DOUBLE)
        _bse_singlet_coefficients.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT + AmB) * eigenvectors.col(level));
        _bse_singlet_coefficients_AR.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT - AmB) * eigenvectors.col(level));
#else
        _bse_singlet_coefficients.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT + AmB) * eigenvectors.col(level)).cast<float>();
        _bse_singlet_coefficients_AR.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT - AmB) * eigenvectors.col(level)).cast<float>();
#endif

      }

      return;
    }
        
    void BSE::printFragInfo(const Population& pop, int i){
      CTP_LOG(ctp::logINFO, _log) << format("           Fragment A -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
              % (100.0 * pop.popH[i](0)) % (100.0 * pop.popE[i](0)) % (pop.Crgs[i](0)) % (pop.Crgs[i](0) + pop.popGs(0)) << flush;
      CTP_LOG(ctp::logINFO, _log) << format("           Fragment B -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
              % (100.0 * pop.popH[i](1)) % (100.0 * pop.popE[i](1)) % (pop.Crgs[i](1)) % (pop.Crgs[i](1) + pop.popGs(1)) << flush;
      return;
    }

    void BSE::printWeights(int i_bse, double weight){
        
      vc2index vc=vc2index(_opt.vmin,_bse_cmin,_bse_ctotal);
      if (weight > _opt.min_print_weight) {
        CTP_LOG(ctp::logINFO, _log) << format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%")
                % (_opt.homo - vc.v(i_bse)) % (vc.c(i_bse) - _opt.homo - 1) % (100.0 * weight) << flush;
      }
      return;
    }

    void BSE::Analyze_singlets(const AOBasis& dftbasis) {

      Interaction act;
      Population pop;
      QMStateType singlet=QMStateType(QMStateType::Singlet);
      std::vector< tools::vec > transition_dipoles=CalcCoupledTransition_Dipoles(dftbasis);
      _orbitals.TransitionDipoles()=transition_dipoles;
      std::vector<double> oscs = _orbitals.Oscillatorstrengths();
      
      if(tools::globals::verbose){
        act = Analyze_eh_interaction(singlet);
      }
      if(dftbasis.getAOBasisFragA() > 0 && dftbasis.getAOBasisFragB()>0){
        pop=FragmentPopulations(singlet,dftbasis);
        _orbitals.setFragmentChargesSingEXC(pop.Crgs);
        _orbitals.setFragment_E_localisation_singlet(pop.popE);
        _orbitals.setFragment_H_localisation_singlet(pop.popH);
        _orbitals.setFragmentChargesGS(pop.popGs);
      }
      
      double hrt2ev = tools::conv::hrt2ev;
      CTP_LOG(ctp::logINFO, _log) << "  ====== singlet energies (eV) ====== "<< flush;
      for (int i = 0; i < _opt.nmax; ++i) {
        const tools::vec& trdip = transition_dipoles[i];
        double osc = oscs[i];
        if (tools::globals::verbose) {
          CTP_LOG(ctp::logINFO, _log) << format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_x> = %5$+1.4f <K_d> = %6$+1.4f")
                  % (i + 1) % (hrt2ev * _bse_singlet_energies(i)) % (1240.0 / (hrt2ev * _bse_singlet_energies(i)))
                  % (hrt2ev * act.qp_contrib(i)) % (hrt2ev * act.exchange_contrib(i)) % (hrt2ev * act.direct_contrib(i)) << flush;
        } else {
          CTP_LOG(ctp::logINFO, _log) << format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                  % (i + 1) % (hrt2ev * _bse_singlet_energies(i)) % (1240.0 / (hrt2ev * _bse_singlet_energies(i))) << flush;
        }
        CTP_LOG(ctp::logINFO, _log) << format("           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f")
                % trdip.getX() % trdip.getY() % trdip.getZ() % (trdip * trdip) % osc << flush;
        for (int i_bse = 0; i_bse < _bse_size; ++i_bse) {
          // if contribution is larger than 0.2, print
          double weight = std::pow(_bse_singlet_coefficients(i_bse, i), 2);
          if (_bse_singlet_coefficients_AR.rows()>0){
               weight-= std::pow(_bse_singlet_coefficients_AR(i_bse, i), 2);
          }
          printWeights(i_bse, weight);
        }
        // results of fragment population analysis 
        if (dftbasis.getAOBasisFragA() > 0 && dftbasis.getAOBasisFragB()>0) {
          printFragInfo(pop, i);
        }

        CTP_LOG(ctp::logINFO, _log) << flush;
      }
      return;
    }
    
    void BSE::Analyze_triplets(const AOBasis& dftbasis) {

      Interaction act;
      Population pop;
      QMStateType triplet=QMStateType(QMStateType::Triplet);
      if(tools::globals::verbose){
        act = Analyze_eh_interaction(triplet);
      }
      if(dftbasis.getAOBasisFragA()>0){
        pop=FragmentPopulations(triplet,dftbasis);
        _orbitals.setFragmentChargesTripEXC(pop.Crgs);
        _orbitals.setFragment_E_localisation_triplet(pop.popE);
        _orbitals.setFragment_H_localisation_triplet(pop.popH);
        _orbitals.setFragmentChargesGS(pop.popGs);
      }
      CTP_LOG(ctp::logINFO, _log) << "  ====== triplet energies (eV) ====== " << flush;
      for (int i = 0; i < _opt.nmax; ++i) {
        if (tools::globals::verbose) {
          CTP_LOG(ctp::logINFO, _log) << format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_d> = %5$+1.4f")
                  % (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) % (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i)))
                  % (tools::conv::hrt2ev * act.qp_contrib(i)) % (tools::conv::hrt2ev *act.direct_contrib(i)) << flush;
        } else {
          CTP_LOG(ctp::logINFO, _log) << format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                  % (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) % (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i))) << flush;
        }
        for (int i_bse = 0; i_bse < _bse_size; ++i_bse) {
          // if contribution is larger than 0.2, print
          double weight = pow(_bse_triplet_coefficients(i_bse, i), 2);
          printWeights(i_bse, weight);
        }
        // results of fragment population analysis 
        if (dftbasis.getAOBasisFragA() > 0) {
          printFragInfo(pop, i);
        }
        CTP_LOG(ctp::logINFO, _log) << format("   ") << flush;
      }
      // storage to orbitals object

      return;
    }

    Eigen::VectorXd BSE::Analyze_IndividualContribution(const QMStateType& type, const BSE_OPERATOR& H){
        Eigen::VectorXd contrib=Eigen::VectorXd::Zero(_opt.nmax);
        if (type == QMStateType::Singlet) {
            for (int i_exc = 0; i_exc < _opt.nmax; i_exc++) {
                MatrixXfd slice_R = _bse_singlet_coefficients.block(0, i_exc, _bse_size, 1);
                contrib(i_exc) =  (slice_R.transpose()*(H * slice_R)).value();
                if (_bse_singlet_coefficients_AR.cols() > 0) {
                    MatrixXfd slice_AR = _bse_singlet_coefficients_AR.block(0, i_exc, _bse_size, 1);
                    // get anti-resonant contribution from direct Keh 
                    contrib(i_exc)-= (slice_AR.transpose()*(H * slice_AR)).value();           
                }
            }
        } else if (type == QMStateType::Triplet) {
            for (int i_exc = 0; i_exc < _opt.nmax; i_exc++) {
                MatrixXfd _slice_R = _bse_triplet_coefficients.block(0, i_exc, _bse_size, 1);
                contrib(i_exc) =  (_slice_R.transpose()*(H * _slice_R)).value();
            }
        } else {
            throw std::runtime_error("BSE::Analyze_eh_interaction:Spin not known!");
        }
        return contrib;
    }


    BSE::Interaction BSE::Analyze_eh_interaction(const QMStateType& type) {

      Interaction analysis;

      BSE_OPERATOR h (_orbitals, _log, _Mmn, _Hqp);
      h.setHqp(1);
      //Eigen::MatrixXd H = h.get_full_matrix();
      analysis.qp_contrib=Analyze_IndividualContribution(type,h);
      
      h.setHqp(0);
      h.setHd(1);
      //H = h.get_full_matrix();
      analysis.direct_contrib=Analyze_IndividualContribution(type,h);

      if (type == QMStateType::Singlet) {
          h.setHd(0);
          h.setHx(1);
          //H = h.get_full_matrix();
          analysis.exchange_contrib=Analyze_IndividualContribution(type,h);
      } else {
            analysis.exchange_contrib=Eigen::VectorXd::Zero(0);
      }
      
      return analysis;
    }

    BSE::Population BSE::FragmentPopulations(const QMStateType& type, const AOBasis& dftbasis) {
      Population pop;
      // Mulliken fragment population analysis
        AOOverlap dftoverlap;
        dftoverlap.Fill(dftbasis);
        CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << dftoverlap.Matrix().rows() << flush;
        // ground state populations
        Eigen::MatrixXd DMAT = _orbitals.DensityMatrixGroundState();
        Eigen::VectorXd nuccharges = _orbitals.FragmentNuclearCharges(dftbasis.getAOBasisFragA());
        Eigen::VectorXd pops = _orbitals.LoewdinPopulation(DMAT, dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
        pop.popGs=nuccharges - pops;
        // population to electron charges and add nuclear charges         
        for (int i_state = 0; i_state < _opt.nmax; i_state++) {
          QMState state=QMState(type,i_state,false);
          // checking Density Matrices
          std::vector< Eigen::MatrixXd > DMAT = _orbitals.DensityMatrixExcitedState(state);
          // hole part
          Eigen::VectorXd popsH = _orbitals.LoewdinPopulation(DMAT[0], dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
          pop.popH.push_back(popsH);
          // electron part
          Eigen::VectorXd popsE = _orbitals.LoewdinPopulation(DMAT[1], dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
          pop.popE.push_back(popsE);
          // update effective charges
          Eigen::VectorXd diff = popsH - popsE;
          pop.Crgs.push_back(diff);
        }
        CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Ran Excitation fragment population analysis " << flush;
     
      return pop;
    }

    std::vector<Eigen::MatrixXd > BSE::CalcFreeTransition_Dipoles(const AOBasis& dftbasis) {
      const Eigen::MatrixXd& dft_orbitals = _orbitals.MOCoefficients();
      // Testing electric dipole AOMatrix
      AODipole dft_dipole;
      dft_dipole.Fill(dftbasis);

      // now transition dipole elements for free interlevel transitions
      std::vector<Eigen::MatrixXd > interlevel_dipoles;

      Eigen::MatrixXd empty = dft_orbitals.block(0,_bse_cmin,dftbasis.AOBasisSize() , _bse_ctotal);
      Eigen::MatrixXd occ = dft_orbitals.block(0,_opt.vmin, dftbasis.AOBasisSize(), _bse_vtotal);
      for (int i_comp = 0; i_comp < 3; i_comp++) {
        interlevel_dipoles.push_back(occ.transpose() * dft_dipole.Matrix()[i_comp] * empty);
      }
      return interlevel_dipoles;
    }

    std::vector<tools::vec > BSE::CalcCoupledTransition_Dipoles(const AOBasis& dftbasis) {
    std::vector<Eigen::MatrixXd > interlevel_dipoles= CalcFreeTransition_Dipoles(dftbasis);
    vc2index vc=vc2index(0,0,_bse_ctotal);
    std::vector<tools::vec > dipols;
    const double sqrt2 = sqrt(2.0);
      for (int i_exc = 0; i_exc < _opt.nmax; i_exc++) {
        tools::vec tdipole = tools::vec(0, 0, 0);
        for (int c = 0; c < _bse_ctotal; c++) {
          for (int v = 0; v < _bse_vtotal; v++) {
            int index_vc = vc.I(v,c);
            double factor = _bse_singlet_coefficients(index_vc, i_exc);
            if (_bse_singlet_coefficients_AR.rows()>0) {
              factor += _bse_singlet_coefficients_AR(index_vc, i_exc);
            }
            // The Transition dipole is sqrt2 bigger because of the spin, the excited state is a linear combination of 2 slater determinants, where either alpha or beta spin electron is excited
            tdipole.x() += factor * interlevel_dipoles[0](v, c);
            tdipole.y() += factor * interlevel_dipoles[1](v, c);
            tdipole.z() += factor * interlevel_dipoles[2](v, c);
          }
        }
        
        dipols.push_back(-sqrt2 *tdipole);     //- because electrons are negative
      }
      return dipols;
    }
    
  }
};
