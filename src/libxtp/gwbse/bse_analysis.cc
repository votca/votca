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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/gwbse.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // MBPT MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //
        
        void GWBSE::BSE_analyze_singlets_BTDA(AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals) {

            
            // expectation values, contributions from e-h coupling
            std::vector<real_gwbse> _contrib_x(_bse_nprint, 0.0);
            std::vector<real_gwbse> _contrib_d(_bse_nprint, 0.0);
            std::vector<real_gwbse> _contrib_qp(_bse_nprint, 0.0);
            BSE_analyze_electron_hole_interaction_BTDA( _bse_singlet_energies, _bse_singlet_coefficients, _bse_singlet_coefficients_AR, _contrib_x, _contrib_d, _contrib_qp );
             
            /*
            // fragment populations
            std::vector<double> _popHA;
            std::vector<double> _popHB;
            std::vector<double> _popEA;
            std::vector<double> _popEB;
              // excitation populations
            std::vector<double> &_CrgsA = _orbitals->FragmentAChargesSingEXC();
            std::vector<double> &_CrgsB = _orbitals->FragmentBChargesSingEXC();
            BSE_FragmentPopulations(dftbasis, _dft_orbitals, _bse_singlet_coefficients, _popHA, _popEA, _popHB, _popEA, _CrgsA, _CrgsB);
            */
            
            
            // for transition dipole moments
            BSE_FreeTransition_Dipoles(dftbasis, _dft_orbitals);
            std::vector<double> _oscillator_strength;
            std::vector<double> _transition_dipole_strength;
            BSE_CoupledTransition_Dipoles_BTDA(_bse_singlet_energies, _bse_singlet_coefficients, _bse_singlet_coefficients_AR, _oscillator_strength, _transition_dipole_strength);
           

            // REPORTING
            // read ground state fragment charges from orbitals object
            //double &_popA = _orbitals->FragmentAChargesGS();
            //double &_popB = _orbitals->FragmentBChargesGS();
            std::vector<ub::vector<double> >& _transition_dipoles = _orbitals->TransitionDipoles();
            
            LOG(ctp::logINFO, *_pLog) << (format("  ====== singlet energies (eV) ====== ")).str() << flush;
            for (int _i = 0; _i < _bse_nprint; _i++) {
                if (votca::tools::globals::verbose) {
                    LOG(ctp::logINFO, *_pLog) << (format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_x> = %5$+1.4f <K_d> = %6$+1.4f")
                            % (_i + 1) % (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)) % (1240.0 / (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)))
                            % (votca::tools::conv::ryd2ev * _contrib_qp[_i]) % (votca::tools::conv::ryd2ev * _contrib_x[_i]) % (votca::tools::conv::ryd2ev * _contrib_d[ _i ])).str() << flush;
                } else {
                    LOG(ctp::logINFO, *_pLog) << (format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                            % (_i + 1) % (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)) % (1240.0 / (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)))).str() << flush;
                }
                
                LOG(ctp::logINFO, *_pLog) << (format("           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f")
                        % (_transition_dipoles[_i](0)) % (_transition_dipoles[_i](1)) % (_transition_dipoles[_i](2)) % (_transition_dipole_strength[_i])
                        % (_oscillator_strength[_i])).str() << flush;
                for (unsigned _i_bse = 0; _i_bse < _bse_size; _i_bse++) {
                    // if contribution is larger than 0.2, print
                    double _weight = pow(_bse_singlet_coefficients(_i_bse, _i), 2) -  pow(_bse_singlet_coefficients_AR(_i_bse, _i), 2);
                    if (_weight > 0.2) {
                        LOG(ctp::logINFO, *_pLog) << (format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%")
                                % (_homo - _index2v[_i_bse]) % (_index2c[_i_bse] - _homo - 1) % (100.0 * _weight)).str() << flush;
                    }
                }
                 
                // results of fragment population analysis 
                /* if (_fragA > 0) {
                    LOG(ctp::logINFO, *_pLog) << (format("           Fragment A -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
                            % (100.0 * _popHA[_i]) % (100.0 * _popEA[_i]) % (_CrgsA[_i]) % (_CrgsA[_i] + _popA)).str() << flush;
                    LOG(ctp::logINFO, *_pLog) << (format("           Fragment B -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
                            % (100.0 * _popHB[_i]) % (100.0 * _popEB[_i]) % (_CrgsB[_i]) % (_CrgsB[_i] + _popB)).str() << flush;
                }*/
                LOG(ctp::logINFO, *_pLog) << (format("   ")).str() << flush;
            }

            // storage to orbitals object
            if (_store_bse_singlets) {

                _bse_singlet_energies.resize(_bse_nmax);
                _bse_singlet_coefficients.resize(_bse_size, _bse_nmax);
                _bse_singlet_coefficients_AR.resize(_bse_size, _bse_nmax);
                //_transition_dipoles.resize(_bse_nprint);
            } else {
                _bse_singlet_energies.resize(0);
                _bse_singlet_coefficients.resize(0, 0);
                _bse_singlet_coefficients_AR.resize(0, 0);
                //_transition_dipoles.resize(0);
            }


        } 

        
        
        
         void GWBSE::BSE_analyze_electron_hole_interaction_BTDA(ub::vector<real_gwbse>& _bse_energies, ub::matrix<real_gwbse>& _bse_coefficients, ub::matrix<real_gwbse>& _bse_coefficients_AR, std::vector<real_gwbse>& _c_x, std::vector<real_gwbse>& _c_d, std::vector<real_gwbse>& _c_qp) {


            if (votca::tools::globals::verbose) {
                for (int _i_exc = 0; _i_exc < _bse_nprint; _i_exc++) {

                    ub::matrix<real_gwbse> _slice_R = ub::project(_bse_coefficients, ub::range(0, _bse_size), ub::range(_i_exc, _i_exc + 1));
                    ub::matrix<real_gwbse> _slice_AR = ub::project(_bse_coefficients_AR, ub::range(0, _bse_size), ub::range(_i_exc, _i_exc + 1));

                    ub::matrix<real_gwbse> _H_R = _eh_d + 2.0 * _eh_x;
                    //ub::matrix<real_gwbse> _H_AR = _eh_d2 + 2.0 * _eh_x;
                    
                    ub::matrix<real_gwbse> tmp = ub::prod( _H_R, _slice_R );
                    ub::matrix<real_gwbse> XAX = ub::prod( ub::trans(_slice_R), tmp  );
                    
                   // tmp = ub::prod( _H_AR, _slice_AR );
                   // ub::matrix<real_gwbse> XBY = ub::prod( ub::trans(_slice_R), tmp  );
                    
                   // tmp = ub::prod( _H_AR, _slice_R );
                   // ub::matrix<real_gwbse> YBX = ub::prod( ub::trans(_slice_AR), tmp  );
                    
                    tmp = ub::prod( _H_R, _slice_AR );
                    ub::matrix<real_gwbse> YAY = ub::prod( ub::trans(_slice_AR), tmp  );
                    
                    /* cout << " Individual "<< endl;
                    cout << "Explicit: " << _bse_energies(_i_exc) << endl;
                    cout << " XAX    : " << XAX(0,0) << endl;
                    //cout << " XBY    : " << XBY(0,0) << endl;
                    //cout << " YBX    : " << YBX(0,0) << endl;
                    cout << " YAY    : " << YAY(0,0) << endl;
                    cout << " SUM    : " << XAX(0,0)+XBY(0,0)-YBX(0,0)-YAY(0,0) << endl;
                    
                    exit(0); */
                    
                    // get resonant contribution from direct Keh 
                    ub::matrix<real_gwbse> _temp = ub::prod(_eh_d - _eh_qp, _slice_R);
                    ub::matrix<real_gwbse> _res  = ub::prod(ub::trans(_slice_R), _temp);
                    _c_d[_i_exc] = _res(0, 0);
                    
                    // get anti-resonant contribution from direct Keh 
                    _temp = ub::prod(_eh_d - _eh_qp, _slice_AR);
                    _res  = ub::prod(ub::trans(_slice_AR), _temp);
                    _c_d[_i_exc] -= _res(0, 0);
                    _c_qp[_i_exc] = _bse_energies(_i_exc) - _c_d[_i_exc] ;
                    
                    

                    // get exchange contribution only, if _c_x is of normal size (singlet)
                    if (_c_x.size() > 0) {
  
                        // resonant part
                        _temp = ub::prod(_eh_x, _slice_R);
                        _res = ub::prod(ub::trans(_slice_R), _temp);
                        _c_x[_i_exc]  = 2.0 * _res(0, 0);
                        
                        // anti-resonant part
                        _temp = ub::prod(_eh_x, _slice_AR);
                        _res = ub::prod(ub::trans(_slice_AR), _temp);
                        _c_x[_i_exc]  -= 2.0 * _res(0, 0);
                        _c_qp[_i_exc] -= _c_x[_i_exc];
                    }
                    
                }
            }

         }
        

        void GWBSE::BSE_analyze_electron_hole_interaction(ub::vector<real_gwbse>& _bse_energies, ub::matrix<real_gwbse>& _bse_coefficients, std::vector<real_gwbse>& _c_x, std::vector<real_gwbse>& _c_d, std::vector<real_gwbse>& _c_qp) {


            if (votca::tools::globals::verbose) {
                for (int _i_exc = 0; _i_exc < _bse_nprint; _i_exc++) {

                    ub::matrix<real_gwbse> _slice = ub::project(_bse_coefficients, ub::range(0, _bse_size), ub::range(_i_exc, _i_exc + 1));

                    ub::matrix<real_gwbse> _temp = ub::prod(_eh_d - _eh_qp, _slice);
                    ub::matrix<real_gwbse> _res  = ub::prod(ub::trans(_slice), _temp);
                    _c_d[_i_exc] = _res(0, 0);
                    _c_qp[_i_exc] = _bse_energies(_i_exc) - _c_d[_i_exc] ;

                    // get exchange contribution only, if _c_x is of normal size (singlet)
                    if (_c_x.size() > 0) {
                        _temp = ub::prod(_eh_x, _slice);
                        _res = ub::prod(ub::trans(_slice), _temp);
                        _c_x[_i_exc]  = 2.0 * _res(0, 0);
                        _c_qp[_i_exc] -= _c_x[_i_exc];
                    }

                }
            }



        }
        void GWBSE::BSE_FragmentPopulations(AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals, ub::matrix<real_gwbse>& _bse_coefficients, std::vector<double>& popHA, std::vector<double>& popEA, std::vector<double>& popHB, std::vector<double>& popEB, std::vector<double> &_CrgsA, std::vector<double> &_CrgsB) {
            double &_popA = _orbitals->FragmentAChargesGS();
            double &_popB = _orbitals->FragmentBChargesGS();

          

            // Mulliken fragment population analysis
            if (_fragA > 0) {

                // get overlap matrix for DFT basisset
                AOOverlap _dftoverlap;
                // initialize overlap matrix
                _dftoverlap.Initialize(dftbasis._AOBasisSize);
                // Fill overlap
                _dftoverlap.Fill(&dftbasis);
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftoverlap._aomatrix.size1() << flush;
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Running Excitation fragment population analysis " << flush;
                // ground state populations
                ub::matrix<double> DMAT = _orbitals->DensityMatrixGroundState(_dft_orbitals);
                double nucA;
                double nucB;
                _orbitals->FragmentNuclearCharges(_fragA, nucA, nucB);
                _orbitals->MullikenPopulation(DMAT, _dftoverlap._aomatrix, dftbasis._AOBasisFragA, _popA, _popB);
                // population to electron charges and add nuclear charges
                _popA = nucA - _popA;
                _popB = nucB - _popB;

                for (int _i_state = 0; _i_state < _bse_nprint; _i_state++) {

                    // checking Density Matrices
                    std::vector<ub::matrix<double> > DMAT = _orbitals->DensityMatrixExcitedState(_dft_orbitals, _bse_coefficients, _i_state);

                    double _totalA;
                    double _totalB;

                    // hole part
                    _orbitals->MullikenPopulation(DMAT[0], _dftoverlap._aomatrix, dftbasis._AOBasisFragA, _totalA, _totalB);

                    popHA.push_back(std::abs(_totalA));
                    popHB.push_back(std::abs(_totalB));

                    // electron part
                    _orbitals->MullikenPopulation(DMAT[1], _dftoverlap._aomatrix, dftbasis._AOBasisFragA, _totalA, _totalB);
                    popEA.push_back(_totalA);
                    popEB.push_back(_totalB);

                    // update effective charges
                    _CrgsA.push_back(popHA[_i_state] - popEA[_i_state]);
                    _CrgsB.push_back(popHB[_i_state] - popEB[_i_state]);

                }

                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " --- complete! " << flush;

            }



        }
        
        void GWBSE::BSE_FreeTransition_Dipoles(AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals){
            
             // Testing electric dipole AOMatrix
            AODipole _dft_dipole;
            _dft_dipole.Initialize(dftbasis._AOBasisSize);
            _dft_dipole.Fill(&dftbasis);
            ub::matrix<double> _temp;
            // now transition dipole elements for free interlevel transitions
            _interlevel_dipoles.resize(3);
            for (int _i_comp = 0; _i_comp < 3; _i_comp++) {
                // 1st: multiply each _dft_momentum component with range of empty states considered in BSE
                ub::matrix<double> _levels = ub::project(_dft_orbitals, ub::range(_bse_cmin, _bse_cmax + 1), ub::range(0, dftbasis._AOBasisSize));
                _temp = ub::prod(_dft_dipole._aomatrix[_i_comp], ub::trans(_levels));
                // 2nd: multiply this with range of occupied states considered in BSE
                _levels = ub::project(_dft_orbitals, ub::range(_bse_vmin, _bse_vmax + 1), ub::range(0, dftbasis._AOBasisSize));
                _interlevel_dipoles[ _i_comp ] = ub::prod(_levels, _temp);
            }
            _temp.resize(0, 0);
            _dft_dipole.Cleanup();
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculated free interlevel transition dipole moments " << flush;
            
            
            
        }

        void GWBSE::BSE_CoupledTransition_Dipoles_BTDA(ub::vector<real_gwbse>& _bse_energies,
                ub::matrix<real_gwbse>& _bse_coefficients,
                ub::matrix<real_gwbse>& _bse_coefficients_AR,
                std::vector<double>& _oscillator_strength,
                std::vector<double>& _transition_dipole_strength) {

            double sqrt2 = sqrt(2.0);
            std::vector<ub::vector<double> >& _transition_dipoles = _orbitals->TransitionDipoles();
            for (int _i_exc = 0; _i_exc < _bse_nprint; _i_exc++) {
                ub::vector<double> _tdipole = ub::zero_vector<double>(3);

                for (unsigned _v = 0; _v < _bse_vtotal; _v++) {
                    for (unsigned _c = 0; _c < _bse_ctotal; _c++) {

                        int index_vc = _bse_ctotal * _v + _c;

                        for (int _i_comp = 0; _i_comp < 3; _i_comp++) {
                            // The Transition dipole is sqrt2 bigger because of the spin, the excited state is a linear combination of 2 slater determinants, where either alpha or beta spin electron is excited
                            _tdipole(_i_comp) += sqrt2 * _bse_singlet_coefficients(index_vc, _i_exc) * _interlevel_dipoles[_i_comp](_v, _c);
                            // anti-resonant part
                            _tdipole(_i_comp) += sqrt2 * _bse_singlet_coefficients_AR(index_vc, _i_exc) * _interlevel_dipoles[_i_comp](_v, _c);
                        }

                    }

                }
                //                cout << _i_exc << " : " << _tdipole[0] << " : "  << _tdipole[1] << " : " << _tdipole[2] << endl;
                _transition_dipoles.push_back(_tdipole);
                _transition_dipole_strength.push_back((ub::inner_prod(_tdipole, _tdipole)));
                //cout<< ub::inner_prod(_tdipole,_tdipole)<< " "<< _tdipole<< endl;
                _oscillator_strength.push_back(_transition_dipole_strength[_i_exc] / 3.0 * (_bse_singlet_energies(_i_exc)));
            }


        }
        
        
        
        
        void GWBSE::BSE_CoupledTransition_Dipoles(ub::vector<real_gwbse>& _bse_energies,
                ub::matrix<real_gwbse>& _bse_coefficients,
                std::vector<double>& _oscillator_strength,
                std::vector<double>& _transition_dipole_strength) {

            double sqrt2 = sqrt(2.0);
            std::vector<ub::vector<double> >& _transition_dipoles = _orbitals->TransitionDipoles();
            for (int _i_exc = 0; _i_exc < _bse_nprint; _i_exc++) {
                ub::vector<double> _tdipole = ub::zero_vector<double>(3);

                for (unsigned _v = 0; _v < _bse_vtotal; _v++) {
                    for (unsigned _c = 0; _c < _bse_ctotal; _c++) {

                        int index_vc = _bse_ctotal * _v + _c;

                        for (int _i_comp = 0; _i_comp < 3; _i_comp++) {
                            // The Transition dipole is sqrt2 bigger because of the spin, the excited state is a linear combination of 2 slater determinants, where either alpha or beta spin electron is excited
                            _tdipole(_i_comp) += sqrt2 * _bse_singlet_coefficients(index_vc, _i_exc) * _interlevel_dipoles[_i_comp](_v, _c);
                        }

                    }

                }
                //                cout << _i_exc << " : " << _tdipole[0] << " : "  << _tdipole[1] << " : " << _tdipole[2] << endl;
                _transition_dipoles.push_back(_tdipole);
                _transition_dipole_strength.push_back((ub::inner_prod(_tdipole, _tdipole)));
                //cout<< ub::inner_prod(_tdipole,_tdipole)<< " "<< _tdipole<< endl;
                _oscillator_strength.push_back(_transition_dipole_strength[_i_exc] / 3.0 * (_bse_singlet_energies(_i_exc)));
            }


        }
        
        void GWBSE::BSE_analyze_singlets(AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals) {


            // expectation values, contributions from e-h coupling
            std::vector<real_gwbse> _contrib_x(_bse_nprint, 0.0);
            std::vector<real_gwbse> _contrib_d(_bse_nprint, 0.0);
            std::vector<real_gwbse> _contrib_qp(_bse_nprint, 0.0);
            BSE_analyze_electron_hole_interaction( _bse_singlet_energies, _bse_singlet_coefficients, _contrib_x, _contrib_d, _contrib_qp );

            // fragment populations
            std::vector<double> _popHA;
            std::vector<double> _popHB;
            std::vector<double> _popEA;
            std::vector<double> _popEB;
              // excitation populations
            std::vector<double> &_CrgsA = _orbitals->FragmentAChargesSingEXC();
            std::vector<double> &_CrgsB = _orbitals->FragmentBChargesSingEXC();
            BSE_FragmentPopulations(dftbasis, _dft_orbitals, _bse_singlet_coefficients, _popHA, _popEA, _popHB, _popEA, _CrgsA, _CrgsB);
            
            // for transition dipole moments
            BSE_FreeTransition_Dipoles(dftbasis, _dft_orbitals);
            std::vector<double> _oscillator_strength;
            std::vector<double> _transition_dipole_strength;
            BSE_CoupledTransition_Dipoles(_bse_singlet_energies, _bse_singlet_coefficients, _oscillator_strength, _transition_dipole_strength);
           

            // REPORTING
            // read ground state fragment charges from orbitals object
            double &_popA = _orbitals->FragmentAChargesGS();
            double &_popB = _orbitals->FragmentBChargesGS();
            std::vector<ub::vector<double> >& _transition_dipoles = _orbitals->TransitionDipoles();
            LOG(ctp::logINFO, *_pLog) << (format("  ====== singlet energies (eV) ====== ")).str() << flush;
            for (int _i = 0; _i < _bse_nprint; _i++) {
                if (votca::tools::globals::verbose) {
                    LOG(ctp::logINFO, *_pLog) << (format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_x> = %5$+1.4f <K_d> = %6$+1.4f")
                            % (_i + 1) % (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)) % (1240.0 / (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)))
                            % (votca::tools::conv::ryd2ev * _contrib_qp[_i]) % (votca::tools::conv::ryd2ev * _contrib_x[_i]) % (votca::tools::conv::ryd2ev * _contrib_d[ _i ])).str() << flush;
                } else {
                    LOG(ctp::logINFO, *_pLog) << (format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                            % (_i + 1) % (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)) % (1240.0 / (votca::tools::conv::ryd2ev * _bse_singlet_energies(_i)))).str() << flush;
                }
                LOG(ctp::logINFO, *_pLog) << (format("           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f")
                        % (_transition_dipoles[_i](0)) % (_transition_dipoles[_i](1)) % (_transition_dipoles[_i](2)) % (_transition_dipole_strength[_i])
                        % (_oscillator_strength[_i])).str() << flush;
                for (unsigned _i_bse = 0; _i_bse < _bse_size; _i_bse++) {
                    // if contribution is larger than 0.2, print
                    double _weight = pow(_bse_singlet_coefficients(_i_bse, _i), 2);
                    if (_weight > 0.2) {
                        LOG(ctp::logINFO, *_pLog) << (format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%")
                                % (_homo - _index2v[_i_bse]) % (_index2c[_i_bse] - _homo - 1) % (100.0 * _weight)).str() << flush;
                    }
                }
                // results of fragment population analysis 
                if (_fragA > 0) {
                    LOG(ctp::logINFO, *_pLog) << (format("           Fragment A -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
                            % (100.0 * _popHA[_i]) % (100.0 * _popEA[_i]) % (_CrgsA[_i]) % (_CrgsA[_i] + _popA)).str() << flush;
                    LOG(ctp::logINFO, *_pLog) << (format("           Fragment B -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f")
                            % (100.0 * _popHB[_i]) % (100.0 * _popEB[_i]) % (_CrgsB[_i]) % (_CrgsB[_i] + _popB)).str() << flush;
                }
                LOG(ctp::logINFO, *_pLog) << (format("   ")).str() << flush;
            }

            // storage to orbitals object
            if (_store_bse_singlets) {

                _bse_singlet_energies.resize(_bse_nmax);
                _bse_singlet_coefficients.resize(_bse_size, _bse_nmax);
                _transition_dipoles.resize(_bse_nprint);
            } else {
                _bse_singlet_energies.resize(0);
                _bse_singlet_coefficients.resize(0, 0);
                _transition_dipoles.resize(0);
            }


        } 

    
        void GWBSE::BSE_analyze_triplets(AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals) {

            // expectation values, contributions from e-h coupling for _n_print lowest excitations
            std::vector<real_gwbse> _contrib_x(0, 0.0);
            std::vector<real_gwbse> _contrib_d(_bse_nprint, 0.0);
            std::vector<real_gwbse> _contrib_qp(_bse_nprint, 0.0);
            BSE_analyze_electron_hole_interaction(_bse_triplet_energies, _bse_triplet_coefficients, _contrib_x, _contrib_d, _contrib_qp);

            // fragment populations
            std::vector<double> _popHA;
            std::vector<double> _popHB;
            std::vector<double> _popEA;
            std::vector<double> _popEB;
            // excitation populations
            std::vector<double> &_CrgsA = _orbitals->FragmentAChargesTripEXC();
            std::vector<double> &_CrgsB = _orbitals->FragmentBChargesTripEXC();
            BSE_FragmentPopulations(dftbasis, _dft_orbitals, _bse_triplet_coefficients, _popHA, _popEA, _popHB, _popEA, _CrgsA, _CrgsB);

            // REPORTING
            double &_popA = _orbitals->FragmentAChargesGS();
            double &_popB = _orbitals->FragmentBChargesGS();
            LOG(ctp::logINFO, *_pLog) << (format("  ====== triplet energies (eV) ====== ")).str() << flush;
            for (int _i = 0; _i < _bse_nprint; _i++) {

                if (votca::tools::globals::verbose) {
                    LOG(ctp::logINFO, *_pLog) << (format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_d> = %6$+1.4f")
                            % (_i + 1) % (votca::tools::conv::ryd2ev * _bse_triplet_energies(_i)) % (1240.0 / (13.6058 * _bse_triplet_energies(_i)))
                            % (votca::tools::conv::ryd2ev * _contrib_qp[_i]) % (votca::tools::conv::ryd2ev * _contrib_d[ _i ])).str() << flush;
                } else {
                    LOG(ctp::logINFO, *_pLog) << (format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                            % (_i + 1) % (votca::tools::conv::ryd2ev * _bse_triplet_energies(_i)) % (1240.0 / (13.6058 * _bse_triplet_energies(_i)))).str() << flush;
                }

                for (unsigned _i_bse = 0; _i_bse < _bse_size; _i_bse++) {
                    // if contribution is larger than 0.2, print
                    real_gwbse _weight = pow(_bse_triplet_coefficients(_i_bse, _i), 2);
                    if (_weight > 0.2) {
                        LOG(ctp::logINFO, *_pLog) << (format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") % (_homo - _index2v[_i_bse]) % (_index2c[_i_bse] - _homo - 1) % (100.0 * _weight)).str() << flush;
                    }
                }
                // results of fragment population analysis 
                if (_fragA > 0) {
                    LOG(ctp::logINFO, *_pLog) << (format("           Fragment A -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f") % (100.0 * _popHA[_i]) % (100.0 * _popEA[_i]) % (_CrgsA[_i]) % (_CrgsA[_i] + _popA)).str() << flush;
                    LOG(ctp::logINFO, *_pLog) << (format("           Fragment B -- hole: %1$5.1f%%  electron: %2$5.1f%%  dQ: %3$+5.2f  Qeff: %4$+5.2f") % (100.0 * _popHB[_i]) % (100.0 * _popEB[_i]) % (_CrgsB[_i]) % (_CrgsB[_i] + _popB)).str() << flush;
                }
                LOG(ctp::logINFO, *_pLog) << (format("   ")).str() << flush;
            }

            // storage to orbitals object
            if (_store_bse_triplets) {
                _bse_triplet_energies.resize(_bse_nmax);
                _bse_triplet_coefficients.resize(_bse_size, _bse_nmax);
            } else {
                _bse_triplet_energies.resize(0);
                _bse_triplet_coefficients.resize(0, 0);
            }

        }


        
   
        
        

    }
    
 
};
