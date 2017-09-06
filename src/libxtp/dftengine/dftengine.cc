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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Overload of uBLAS prod function with MKL/GSL implementations

#include "votca/xtp/aobasis.h"
#include "votca/xtp/qminterface.h"
#include <votca/xtp/dftengine.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/threecenters.h>

#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>
#include <votca/tools/constants.h>

#include <votca/xtp/elements.h>
#include <votca/xtp/diis.h>

#include <votca/ctp/xinteractor.h>
#include <votca/ctp/logger.h>



using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
        // +++++++++++++++++++++++++++++ //
        // DFTENGINE MEMBER FUNCTIONS        //
        // +++++++++++++++++++++++++++++ //

        void DFTENGINE::CleanUp() {

        }

        void DFTENGINE::Initialize(Property* options) {

            string key = Identify();

            // get OpenMP thread number

            _openmp_threads = options->ifExistsReturnElseReturnDefault<int>(key + ".openmp", 0);

            // basis sets
            _dftbasis_name = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".dftbasis");

            if (options->exists(key + ".auxbasis")) {
                _auxbasis_name = options->get(key + ".auxbasis").as<string>();
                _with_RI = true;
            } else {
                _with_RI = false;
            }

            if (options->exists(key + ".ecp")) {
                _ecp_name = options->get(key + ".ecp").as<string>();
                _with_ecp = true;
            } else {
                _with_ecp = false;
            }
            _with_guess = options->ifExistsReturnElseReturnDefault<bool>(key + ".read_guess", false);
            _initial_guess = options->ifExistsReturnElseReturnDefault<string>(key + ".initial_guess", "atom");


            // numerical integrations
            _grid_name = options->ifExistsReturnElseReturnDefault<string>(key + ".integration_grid", "medium");
            _use_small_grid = options->ifExistsReturnElseReturnDefault<bool>(key + ".integration_grid_small", true);
            _grid_name_small = Choosesmallgrid(_grid_name);

            // exchange and correlation as in libXC

            _xc_functional_name = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".xc_functional");

            _numofelectrons = 0;

            if (options->exists(key + ".convergence")) {

                _Econverged = options->ifExistsReturnElseReturnDefault<double>(key + ".convergence.energy", 1e-7);
                _error_converged = options->ifExistsReturnElseReturnDefault<double>(key + ".convergence.error", 1e-7);
                _max_iter = options->ifExistsReturnElseReturnDefault<int>(key + ".convergence.max_iterations", 100);



                if (options->exists(key + ".convergence.method")) {
                    string method = options->get(key + ".convergence.method").as<string>();
                    if (method == "DIIS") {
                        _usediis = true;
                    } else if (method == "mixing") {
                        _usediis = false;
                    } else {
                        cout << "WARNING method not known. Using Mixing" << endl;
                        _usediis = false;
                    }
                } else {
                    _usediis = true;
                }
                if (!_usediis) {
                    _histlength = 1;
                    _maxout = false;
                }

                if (options->exists(key + ".convergence.mixing")) {
                    _useautomaticmixing = false;
                    _mixingparameter = options->get(key + ".convergence.mixing").as<double>();
                } else {
                    _useautomaticmixing = true;
                    _mixingparameter = -10;
                }

                _levelshift = options->ifExistsReturnElseReturnDefault<double>(key + ".convergence.levelshift", 0.0);
                _levelshiftend = options->ifExistsReturnElseReturnDefault<double>(key + ".convergence.levelshift_end", 0.8);
                _maxout = options->ifExistsReturnElseReturnDefault<bool>(key + ".convergence.DIIS_maxout", false);
                _histlength = options->ifExistsReturnElseReturnDefault<int>(key + ".convergence.DIIS_length", 10);
                _diis_start = options->ifExistsReturnElseReturnDefault<double>(key + ".convergence.DIIS_start", 0.01);
                _adiis_start = options->ifExistsReturnElseReturnDefault<double>(key + ".convergence.ADIIS_start", 2);

            } else {
                _Econverged = 1e-7;
                _error_converged = 1e-7;
                _maxout = false;
                _diis_start = 0.01;
                _adiis_start = 2;
                _histlength = 10;
                _useautomaticmixing = true;
                _mixingparameter = -10;
                _usediis = true;
                _max_iter = 100;
                _levelshift = 0.25;
                _levelshiftend = 0.8;
            }

            return;
        }

        /*
         *    Density Functional theory implementation
         *
         */




        bool DFTENGINE::Evaluate(Orbitals* _orbitals) {


            // set the parallelization
#ifdef _OPENMP

            omp_set_num_threads(_openmp_threads);

#endif





            /**** END OF PREPARATION ****/

            /**** Density-independent matrices ****/

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Inside Evaluate " << flush;


            ub::vector<double>& MOEnergies = _orbitals->MOEnergies();
            ub::matrix<double>& MOCoeff = _orbitals->MOCoefficients();
            if (MOEnergies.size() != _dftbasis.AOBasisSize()) {
                MOEnergies.resize(_dftbasis.AOBasisSize());
            }
            if (MOCoeff.size1() != _dftbasis.AOBasisSize() || MOCoeff.size2() != _dftbasis.AOBasisSize()) {
                MOCoeff.resize(_dftbasis.AOBasisSize(), _dftbasis.AOBasisSize());
            }

            /**** Construct initial density  ****/

            ub::matrix<double> H0 = _dftAOkinetic.Matrix() + _dftAOESP.getNuclearpotential();

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Constructed inital density " << flush;


            NuclearRepulsion();
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Nuclear Rep " << flush;


            if (_addexternalsites) {
                H0 += _dftAOESP.getExternalpotential();

                H0 += _dftAODipole_Potential.getExternalpotential();
                H0 += _dftAOQuadrupole_Potential.getExternalpotential();

                double estat = ExternalRepulsion();
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " E_electrostatic " << estat << flush;
                E_nucnuc += estat;

            }

            if (_do_externalfield) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Integrated external potential on grid " << flush;
                double extneralgrid_nucint = ExternalGridRepulsion(_externalgrid_nuc);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Nuclei external potential interaction " << extneralgrid_nucint << " Hartree" << flush;
                H0 += _gridIntegration_ext.IntegrateExternalPotential(_externalgrid);
                E_nucnuc += extneralgrid_nucint;
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Nuclear Repulsion Energy is " << E_nucnuc << flush;


            // if we have a guess we do not need this.
            if (_with_guess) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Reading guess from orbitals object/file" << flush;
                _dftAOdmat = _orbitals->DensityMatrixGroundState();
            } else if (guess_set) {
                ConfigOrbfile(_orbitals);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using starting guess from last iteration" << flush;
                _dftAOdmat = last_dmat;
            } else {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup Initial Guess using: " << _initial_guess << flush;
                // this temp is necessary because eigenvalues_general returns MO^T and not MO
                if (_initial_guess == "independent") {
                    ub::matrix<double> copy = H0;
                    _diis.SolveFockmatrix(MOEnergies, MOCoeff, copy);
                    _dftAOdmat = _orbitals->DensityMatrixGroundState();

                } else if (_initial_guess == "atom") {

                    _dftAOdmat = AtomicGuess(_orbitals);
                    //cout<<_dftAOdmat<<endl;

                    if (_with_RI) {
                        _ERIs.CalculateERIs(_dftAOdmat);
                    } else {
                        _ERIs.CalculateERIs_4c_small_molecule(_dftAOdmat);
                    }
                    if (_use_small_grid) {
                        _orbitals->AOVxc() = _gridIntegration_small.IntegrateVXC(_dftAOdmat);
                        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled approximate DFT Vxc matrix " << flush;
                    } else {
                        _orbitals->AOVxc() = _gridIntegration.IntegrateVXC(_dftAOdmat);
                        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Vxc matrix " << flush;
                    }
                    ub::matrix<double> H = H0 + _ERIs.getERIs() + _orbitals->AOVxc();
                    _diis.SolveFockmatrix(MOEnergies, MOCoeff, H);
                    _dftAOdmat = _orbitals->DensityMatrixGroundState();
                    //cout<<_dftAOdmat<<endl;
                    //Have to do one full iteration here, levelshift needs MOs;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Full atomic density Matrix gives N=" << std::setprecision(9) << linalg_traceofProd(_dftAOdmat, _dftAOoverlap.Matrix()) << " electrons." << flush;
                } else {
                    throw runtime_error("Initial guess method not known/implemented");
                }

            }


            _orbitals->setQMpackage("xtp");

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " STARTING SCF cycle" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << " --------------------------------------------------------------------------" << flush;

            double energyold = 0;
            double diiserror = 100; //is evolved in DIIs scheme
            Mixing Mixer(_useautomaticmixing, _mixingparameter, &_dftAOoverlap.Matrix(), _pLog);

            for (_this_iter = 0; _this_iter < _max_iter; _this_iter++) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Iteration " << _this_iter + 1 << " of " << _max_iter << flush;

                if (_with_RI) {
                    _ERIs.CalculateERIs(_dftAOdmat);
                } else {
                    _ERIs.CalculateERIs_4c_small_molecule(_dftAOdmat);
                }

                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Electron repulsion matrix of dimension: " << _ERIs.getSize1() << " x " << _ERIs.getSize2() << flush;
                double vxcenergy = 0.0;
                if (_use_small_grid && diiserror > 1e-5) {
                    _orbitals->AOVxc() = _gridIntegration_small.IntegrateVXC(_dftAOdmat);
                    vxcenergy = _gridIntegration_small.getTotEcontribution();
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled approximate DFT Vxc matrix " << flush;
                } else {
                    _orbitals->AOVxc() = _gridIntegration.IntegrateVXC(_dftAOdmat);
                    vxcenergy = _gridIntegration.getTotEcontribution();
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Vxc matrix " << flush;
                }
                //cout<<_dftAOdmat<<endl;
                ub::matrix<double> H = H0 + _ERIs.getERIs() + _orbitals->AOVxc();
                //cout<<H0<<endl;
                //exit(0);
                double Eone = linalg_traceofProd(_dftAOdmat, H0);
                double Etwo = 0.5 * _ERIs.getERIsenergy() + vxcenergy;
                double totenergy = Eone + E_nucnuc + Etwo;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Single particle energy " << std::setprecision(12) << Eone << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Two particle energy " << std::setprecision(12) << Etwo << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << std::setprecision(12) << " Exc contribution " << vxcenergy << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Total Energy " << std::setprecision(12) << totenergy << flush;

                diiserror = _diis.Evolve(_dftAOdmat, H, MOEnergies, MOCoeff, _this_iter, totenergy);
                //cout<<"Energies "<<MOEnergies<<endl;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " DIIs error " << diiserror << flush;

                ub::matrix<double> dmatin = _dftAOdmat;
                _dftAOdmat = _orbitals->DensityMatrixGroundState();
                if (!(diiserror < _adiis_start && _usediis && _this_iter > 2)) {
                    _dftAOdmat = Mixer.MixDmat(dmatin, _dftAOdmat);

                } else {
                    Mixer.Updatemix(dmatin, _dftAOdmat);
                }

                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Updated Density Matrix " << flush;

                if (tools::globals::verbose) {
                    for (int i = 0; i<int(MOEnergies.size()); i++) {
                        if (i <= _numofelectrons / 2 - 1) {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t" << i << " occ " << MOEnergies(i) << flush;
                        } else {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t" << i << " vir " << MOEnergies(i) << flush;

                        }
                    }
                }

                CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\tGAP " << MOEnergies(_numofelectrons / 2) - MOEnergies(_numofelectrons / 2 - 1) << flush;

                if (std::abs(totenergy - energyold) < _Econverged && diiserror < _error_converged) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Total Energy has converged to " << std::setprecision(9) << std::abs(totenergy - energyold) << "[Ha] after " << _this_iter + 1 <<
                            " iterations. DIIS error is converged up to " << _error_converged << "[Ha]" << flush;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Final Single Point Energy " << std::setprecision(12) << totenergy << " Ha" << flush;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " MO Energies  [Ha]" << flush;
                    for (int i = 0; i<int(MOEnergies.size()); i++) {
                        if (i <= _numofelectrons / 2 - 1) {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t" << i << " occ " << std::setprecision(12) << MOEnergies(i) << flush;
                        } else {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t" << i << " vir " << std::setprecision(12) << MOEnergies(i) << flush;
                        }
                    }

                    last_dmat = _dftAOdmat;
                    guess_set = true;
                    // orbitals saves total energies in [eV]
                    _orbitals->setQMEnergy(totenergy * tools::conv::hrt2ev);
                    break;
                } else {
                    energyold = totenergy;
                }
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Density Matrix gives N=" << std::setprecision(9) << linalg_traceofProd(_dftAOdmat, _dftAOoverlap.Matrix()) << " electrons." << flush;

            }
            return true;
        }




        // SETUP INVARIANT AOMatrices

        void DFTENGINE::SetupInvariantMatrices() {


            // local variables for checks
            // check eigenvalues of overlap matrix, if too small basis might have linear dependencies
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;

            {
                // DFT AOOverlap matrix
                _dftAOoverlap.Initialize(_dftbasis.AOBasisSize());
                _dftAOoverlap.Fill(_dftbasis);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftAOoverlap.Dimension() << flush;
                //cout<<"overlap"<<_dftAOoverlap.Matrix()<<endl;
                // check DFT basis for linear dependence
                linalg_eigenvalues(_dftAOoverlap.Matrix(), _eigenvalues, _eigenvectors);

                //Not brilliant but we need S-1/2 for DIIS and I do not want to calculate it each time
                ub::matrix<double> _diagS = ub::zero_matrix<double>(_eigenvectors.size1(), _eigenvectors.size1());
                for (unsigned _i = 0; _i < _eigenvalues.size(); _i++) {

                    _diagS(_i, _i) = 1.0 / sqrt(_eigenvalues[_i]);
                }
                ub::matrix<double> _temp = ub::prod(_diagS, ub::trans(_eigenvectors));
                _Sminusonehalf = ub::prod(_eigenvectors, _temp);

                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Smallest eigenvalue of DFT Overlap matrix : " << _eigenvalues[0] << flush;
            }

            _dftAOkinetic.Initialize(_dftbasis.AOBasisSize());
            _dftAOkinetic.Fill(_dftbasis);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Kinetic energy matrix of dimension: " << _dftAOkinetic.Dimension() << flush;


            _dftAOESP.Initialize(_dftbasis.AOBasisSize());
            _dftAOESP.Fillnucpotential(_dftbasis, _atoms, _with_ecp);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT nuclear potential matrix of dimension: " << _dftAOESP.Dimension() << flush;

            if (_addexternalsites) {
                _dftAOESP.Fillextpotential(_dftbasis, _externalsites);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT external pointcharge potential matrix of dimension: " << _dftAOESP.Dimension() << flush;

                _dftAODipole_Potential.Fillextpotential(_dftbasis, _externalsites);
                if (_dftAODipole_Potential.Dimension() > 0) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT external dipole potential matrix of dimension: " << _dftAODipole_Potential.Dimension() << flush;
                }
                _dftAOQuadrupole_Potential.Fillextpotential(_dftbasis, _externalsites);
                if (_dftAOQuadrupole_Potential.Dimension()) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT external quadrupole potential matrix of dimension: " << _dftAOQuadrupole_Potential.Dimension() << flush;
                }
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " External sites\t Name \t Coordinates \t charge \t dipole \t quadrupole" << flush;


                for (unsigned i = 0; i < _externalsites.size(); i++) {

                    vector<ctp::APolarSite*> ::iterator pit;
                    for (pit = _externalsites[i]->begin(); pit < _externalsites[i]->end(); ++pit) {

                        CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << (*pit)->getName() << " | " << (*pit)->getPos().getX()
                                << " " << (*pit)->getPos().getY() << " " << (*pit)->getPos().getZ() << " | " << (*pit)->getQ00();
                        if ((*pit)->getRank() > 0) {
                            tools::vec dipole = (*pit)->getQ1();
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << " | " << dipole.getX()
                                    << " " << dipole.getY() << " " << dipole.getZ();
                        }
                        if ((*pit)->getRank() > 1) {
                            std::vector<double> quadrupole = (*pit)->getQ2();
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << " | " << quadrupole[0] << " " << quadrupole[1] << " " << quadrupole[2] << " "
                                    << quadrupole[3] << " " << quadrupole[4];
                        }
                        CTP_LOG(ctp::logDEBUG, *_pLog)  << flush;
                    }
                }
            }




            if (_with_ecp) {
                _dftAOECP.Initialize(_dftbasis.AOBasisSize());
                _dftAOECP.Fill(_dftbasis, vec(0, 0, 0), &_ecp);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT ECP matrix of dimension: " << _dftAOECP.Dimension() << flush;
                _dftAOESP.getNuclearpotential() += _dftAOECP.Matrix();
            }

            _diis.Configure(_usediis, true, _histlength, _maxout, _diismethod, _adiis_start, _diis_start, _levelshift, _levelshiftend, _numofelectrons / 2);
            _diis.setLogger(_pLog);
            _diis.setOverlap(&_dftAOoverlap.Matrix());
            _diis.setSqrtOverlap(&_Sminusonehalf);

            if (_with_RI) {

                AOCoulomb _auxAOcoulomb;
                _auxAOcoulomb.Initialize(_auxbasis.AOBasisSize());
                _auxAOcoulomb.Fill(_auxbasis);

                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled AUX Coulomb matrix of dimension: " << _auxAOcoulomb.Dimension() << flush;
                ub::matrix<double> AuxAOcoulomb_inv = ub::zero_matrix<double>(_auxAOcoulomb.Dimension(), _auxAOcoulomb.Dimension());
                int dimensions = linalg_invert_svd(_auxAOcoulomb.Matrix(), AuxAOcoulomb_inv, 1e7);

                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Inverted AUX Coulomb matrix, removed " << dimensions << " functions from aux basis" << flush;


                // prepare invariant part of electron repulsion integrals
                _ERIs.Initialize(_dftbasis, _auxbasis, AuxAOcoulomb_inv);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup invariant parts of Electron Repulsion integrals " << flush;
            } else {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculating 4c integrals. " << flush;
                _ERIs.Initialize_4c_small_molecule(_dftbasis);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculated 4c integrals. " << flush;
            }

            return;
        }

        ub::matrix<double> DFTENGINE::AtomicGuess(Orbitals* _orbitals) {
            ub::matrix<double> guess = ub::zero_matrix<double>(_dftbasis.AOBasisSize());

            std::vector<ctp::QMAtom*> uniqueelements;
            std::vector<ctp::QMAtom*>::const_iterator at;
            std::vector<ctp::QMAtom*>::iterator st;
            std::vector< ub::matrix<double> > uniqueatom_guesses;
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Scanning molecule of size " << _atoms.size() << " for unique elements" << flush;
            for (at = _atoms.begin(); at < _atoms.end(); ++at) {
                bool exists = false;
                if (uniqueelements.size() == 0) {
                    exists = false;
                } else {
                    for (st = uniqueelements.begin(); st < uniqueelements.end(); ++st) {
                        if ((*at)->type == (*st)->type) {
                            exists = true;
                            break;
                        }
                    }
                }
                if (!exists) {
                    uniqueelements.push_back((*at));
                }
            }
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " " << uniqueelements.size() << " unique elements found" << flush;
            Elements _elements;
            for (st = uniqueelements.begin(); st < uniqueelements.end(); ++st) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculating atom density for " << (*st)->type << flush;
                bool with_ecp = _with_ecp;
                if ((*st)->type == "H" || (*st)->type == "He") {
                    with_ecp = false;
                }
                std::vector<ctp::QMAtom*> atom;
                atom.push_back(*st);

                AOBasis dftbasis;
                AOBasis ecp;
                NumericalIntegration gridIntegration;
                dftbasis.AOBasisFill(&_dftbasisset, atom);
                if (with_ecp) {
                    ecp.ECPFill(&_ecpbasisset, atom);
                }
                gridIntegration.GridSetup(_grid_name, &_dftbasisset, atom, &dftbasis);
                gridIntegration.setXCfunctional(_xc_functional_name);
                int numofelectrons = int(_elements.getNucCrg((*st)->type));
                int alpha_e = 0;
                int beta_e = 0;
                if (with_ecp) {
                    numofelectrons -= int(_ecpbasisset.getElement((*st)->type)->getNcore());
                }

                if ((numofelectrons % 2) != 0) {
                    alpha_e = numofelectrons / 2 + numofelectrons % 2;
                    beta_e = numofelectrons / 2;
                } else {
                    alpha_e = numofelectrons / 2;
                    beta_e = alpha_e;
                }

                AOOverlap dftAOoverlap;
                AOKinetic dftAOkinetic;
                AOESP dftAOESP;
                AOECP dftAOECP;
                ub::vector<double> eigenvalues;
                ub::matrix<double> eigenvectors;
                ub::matrix<double> Sminusonehalf;
                ERIs ERIs_atom;

                // DFT AOOverlap matrix
                dftAOoverlap.Initialize(dftbasis.AOBasisSize());
                dftAOoverlap.Fill(dftbasis);
                linalg_eigenvalues(dftAOoverlap.Matrix(), eigenvalues, eigenvectors);

                //Not brilliant but we need S-1/2 for DIIS and I do not want to calculate it each time
                ub::matrix<double> _diagS = ub::zero_matrix<double>(eigenvectors.size1(), eigenvectors.size1());
                for (unsigned _i = 0; _i < eigenvalues.size(); _i++) {

                    _diagS(_i, _i) = 1.0 / sqrt(eigenvalues[_i]);
                }
                ub::matrix<double> _temp = ub::prod(_diagS, ub::trans(eigenvectors));
                Sminusonehalf = ub::prod(eigenvectors, _temp);

                dftAOkinetic.Initialize(dftbasis.AOBasisSize());
                dftAOkinetic.Fill(dftbasis);

                dftAOESP.Initialize(dftbasis.AOBasisSize());
                dftAOESP.Fillnucpotential(dftbasis, atom, with_ecp);
                ERIs_atom.Initialize_4c_small_molecule(dftbasis);

                ub::vector<double>MOEnergies_alpha;
                ub::matrix<double>MOCoeff_alpha;
                ub::vector<double>MOEnergies_beta;
                ub::matrix<double>MOCoeff_beta;
                Diis diis_alpha;
                Diis diis_beta;
                double adiisstart = 0;
                double diisstart = 0.0001;
                Mixing Mix_alpha(true, 0.7, &dftAOoverlap.Matrix(), _pLog);
                Mixing Mix_beta(true, 0.7, &dftAOoverlap.Matrix(), _pLog);
                diis_alpha.Configure(true, false, 20, 0, "", adiisstart, diisstart, 0.1, 0.000, alpha_e);
                diis_alpha.setLogger(_pLog);
                diis_alpha.setOverlap(&dftAOoverlap.Matrix());
                diis_alpha.setSqrtOverlap(&Sminusonehalf);
                diis_beta.Configure(true, false, 20, 0, "", adiisstart, diisstart, 0.1, 0.000, beta_e);
                diis_beta.setLogger(_pLog);
                diis_beta.setOverlap(&dftAOoverlap.Matrix());
                diis_beta.setSqrtOverlap(&Sminusonehalf);
                /**** Construct initial density  ****/

                ub::matrix<double> H0 = dftAOkinetic.Matrix() + dftAOESP.getNuclearpotential();
                if (with_ecp) {
                    dftAOECP.Initialize(dftbasis.AOBasisSize());
                    dftAOECP.Fill(dftbasis, vec(0, 0, 0), &ecp);
                    H0 += dftAOECP.Matrix();
                }
                ub::matrix<double> copy = H0;
                diis_alpha.SolveFockmatrix(MOEnergies_alpha, MOCoeff_alpha, copy);

                MOEnergies_beta = MOEnergies_alpha;
                MOCoeff_beta = MOCoeff_alpha;


                //ub::matrix<double>dftAOdmat_alpha = DensityMatrix_frac(MOCoeff_alpha,MOEnergies_alpha,alpha_e);
                //ub::matrix<double>dftAOdmat_beta = DensityMatrix_frac(MOCoeff_beta,MOEnergies_beta,beta_e);
                ub::matrix<double>dftAOdmat_alpha = DensityMatrix_unres(MOCoeff_alpha, alpha_e);
                if ((*st)->type == "H") {
                    uniqueatom_guesses.push_back(dftAOdmat_alpha);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Atomic density Matrix for " << (*st)->type <<
                            " gives N=" << std::setprecision(9) << linalg_traceofProd(dftAOdmat_alpha, dftAOoverlap.Matrix()) << " electrons." << flush;
                    continue;
                }

                ub::matrix<double>dftAOdmat_beta = DensityMatrix_unres(MOCoeff_beta, beta_e);
                bool _HF = false;
                double energyold = 0;
                int maxiter = 50;
                for (int this_iter = 0; this_iter < maxiter; this_iter++) {

                    ERIs_atom.CalculateERIs_4c_small_molecule(dftAOdmat_alpha + dftAOdmat_beta);
                    double E_two_alpha = linalg_traceofProd(ERIs_atom.getERIs(), dftAOdmat_alpha);
                    double E_two_beta = linalg_traceofProd(ERIs_atom.getERIs(), dftAOdmat_beta);
                    ub::matrix<double> H_alpha = H0 + ERIs_atom.getERIs();
                    ub::matrix<double> H_beta = H0 + ERIs_atom.getERIs();
                    if (_HF) {
                        ERIs_atom.CalculateEXX_4c_small_molecule(dftAOdmat_alpha);
                        double E_exx_alpha = -linalg_traceofProd(ERIs_atom.getEXX(), dftAOdmat_alpha);
                        H_alpha -= ERIs_atom.getEXX();
                        E_two_alpha += E_exx_alpha;
                        ERIs_atom.CalculateEXX_4c_small_molecule(dftAOdmat_beta);
                        double E_exx_beta = -linalg_traceofProd(ERIs_atom.getEXX(), dftAOdmat_beta);
                        H_beta -= ERIs_atom.getEXX();
                        E_two_beta += E_exx_beta;

                    } else {
                        ub::matrix<double> AOVxc_alpha = gridIntegration.IntegrateVXC(dftAOdmat_alpha);
                        double E_vxc_alpha = gridIntegration.getTotEcontribution();
                        ub::matrix<double> AOVxc_beta = gridIntegration.IntegrateVXC(dftAOdmat_beta);
                        double E_vxc_beta = gridIntegration.getTotEcontribution();
                        H_alpha += AOVxc_alpha;
                        H_beta += AOVxc_beta;
                        E_two_alpha += E_vxc_alpha;
                        E_two_beta += E_vxc_beta;

                    }
                    double E_one_alpha = linalg_traceofProd(dftAOdmat_alpha, H0);
                    double E_one_beta = linalg_traceofProd(dftAOdmat_beta, H0);

                    double E_alpha = E_one_alpha + E_two_alpha;
                    double E_beta = E_one_beta + E_two_beta;

                    double totenergy = E_alpha + E_beta;
                    //evolve alpha
                    double diiserror_alpha = diis_alpha.Evolve(dftAOdmat_alpha, H_alpha, MOEnergies_alpha, MOCoeff_alpha, this_iter, E_alpha);

                    ub::matrix<double> dmatin_alpha = dftAOdmat_alpha;
                    dftAOdmat_alpha = DensityMatrix_unres(MOCoeff_alpha, alpha_e);

                    if (!(diiserror_alpha < 0.005 && this_iter > 2)) {
                        dftAOdmat_alpha = Mix_alpha.MixDmat(dmatin_alpha, dftAOdmat_alpha, false);
                        //cout<<"mixing_alpha"<<endl;
                    } else {
                        Mix_alpha.Updatemix(dmatin_alpha, dftAOdmat_alpha);
                    }
                    //evolve beta
                    double diiserror_beta = diis_beta.Evolve(dftAOdmat_beta, H_beta, MOEnergies_beta, MOCoeff_beta, this_iter, E_beta);
                    ub::matrix<double> dmatin_beta = dftAOdmat_beta;
                    dftAOdmat_beta = DensityMatrix_unres(MOCoeff_beta, beta_e);
                    //dftAOdmat_beta=DensityMatrix_frac(MOCoeff_beta,MOEnergies_beta,beta_e);
                    if (!(diiserror_beta < 0.005 && this_iter > 2)) {
                        dftAOdmat_beta = Mix_beta.MixDmat(dmatin_beta, dftAOdmat_beta, false);

                    } else {
                        Mix_beta.Updatemix(dmatin_beta, dftAOdmat_beta);
                    }


                    if (tools::globals::verbose) {
                        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Iter " << this_iter << " of " << maxiter
                                << " Etot " << totenergy << " diise_a " << diiserror_alpha << " diise_b " << diiserror_beta
                                << " a_gap " << MOEnergies_alpha(alpha_e) - MOEnergies_alpha(alpha_e - 1) << " b_gap " << MOEnergies_beta(beta_e) - MOEnergies_beta(beta_e - 1) << flush;
                    }
                    bool converged = (std::abs(totenergy - energyold) < _Econverged && diiserror_alpha < _error_converged && diiserror_beta < _error_converged);
                    if (converged || this_iter == maxiter - 1) {

                        if (converged) {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Converged after " << this_iter + 1 << " iterations" << flush;
                        } else {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Not converged after " << this_iter + 1 <<
                                    " iterations. Using unconverged density. DIIsError_alpha=" << diiserror_alpha << " DIIsError_beta=" << diiserror_beta << flush;
                        }


                        ub::matrix<double> avdmat = AverageShells(dftAOdmat_alpha + dftAOdmat_beta, dftbasis);
                        uniqueatom_guesses.push_back(avdmat);
                        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Atomic density Matrix for " << (*st)->type << " gives N="
                                << std::setprecision(9) << linalg_traceofProd(avdmat, dftAOoverlap.Matrix()) << " electrons." << flush;
                        break;

                    } else {
                        energyold = totenergy;
                    }
                }
            }
            unsigned start = 0;
            unsigned end = 0;
            for (at = _atoms.begin(); at < _atoms.end(); ++at) {
                unsigned index = 0;
                for (unsigned i = 0; i < uniqueelements.size(); i++) {
                    if ((*at)->type == uniqueelements[i]->type) {
                        index = i;
                        break;
                    }
                }

                Element* element = _dftbasisset.getElement((*at)->type);
                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                    end += (*its)->getnumofFunc();
                }
                ub::project(guess, ub::range(start, end), ub::range(start, end)) = uniqueatom_guesses[index];


                start = end;
            }

            return guess;
        }

        void DFTENGINE::ConfigOrbfile(Orbitals* _orbitals) {
            if (_with_guess) {

                if (_orbitals->hasDFTbasis()) {
                    if (_orbitals->getDFTbasis() != _dftbasis_name) {
                        throw runtime_error((boost::format("Basisset Name in guess orb file and in dftengine option file differ %1% vs %2%") % _orbitals->getDFTbasis() % _dftbasis_name).str());
                    }
                } else {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " WARNING: Orbital file has no basisset information,using it as a guess might work or not for calculation with " << _dftbasis_name << flush;
                }
            }
            _orbitals->setDFTbasis(_dftbasis_name);
            _orbitals->setBasisSetSize(_dftbasis.AOBasisSize());
            if (_with_ecp) {
                _orbitals->setECP(_ecp_name);
            }

            if (_with_guess) {
                if(_orbitals->hasECP() || _with_ecp){
                    if(_orbitals->getECP()!=_ecp_name){
                        throw runtime_error((boost::format("ECPs in orb file: %1% and options %2% differ") % _orbitals->getECP() % _ecp_name).str());
                    }
                }
                if (_orbitals->getNumberOfElectrons() != _numofelectrons / 2) {
                    throw runtime_error((boost::format("Number of electron in guess orb file: %1% and in dftengine: %2% differ.") % _orbitals->getNumberOfElectrons() % (_numofelectrons / 2)).str());
                }
                if (_orbitals->getNumberOfLevels() != _dftbasis.AOBasisSize()) {
                    throw runtime_error((boost::format("Number of levels in guess orb file: %1% and in dftengine: %2% differ.") % _orbitals->getNumberOfLevels() % _dftbasis.AOBasisSize()).str());
                }
            } else {
                _orbitals->setNumberOfElectrons(_numofelectrons / 2);
                _orbitals->setNumberOfLevels(_numofelectrons / 2, _dftbasis.AOBasisSize() - _numofelectrons / 2);
            }
            return;
        }


        // PREPARATION

        void DFTENGINE::Prepare(Orbitals* _orbitals) {
            #ifdef _OPENMP

            omp_set_num_threads(_openmp_threads);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using " << omp_get_max_threads() << " threads" << flush;

#endif

            if ( _atoms.size() == 0 ){
            for (const auto& atom : _orbitals->QMAtoms()) {
                if (!atom->from_environment) {
                    _atoms.push_back(atom);
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Molecule Coordinates [A] " << flush;
            for (unsigned i = 0; i < _atoms.size(); i++) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << _atoms[i]->type << " " << _atoms[i]->x << " " << _atoms[i]->y << " " << _atoms[i]->z << " " << flush;
            }

            // load and fill DFT basis set
            _dftbasisset.LoadBasisSet(_dftbasis_name);

            _dftbasis.AOBasisFill(&_dftbasisset, _atoms);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << flush;

            if (_with_RI) {
                // load and fill AUX basis set
                _auxbasisset.LoadBasisSet(_auxbasis_name);
                //_orbitals->setDFTbasis( _dftbasis_name );
                _auxbasis.AOBasisFill(&_auxbasisset, _atoms);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded AUX Basis Set " << _auxbasis_name << flush;
            }
            if (_with_ecp) {
                // load ECP (element-wise information) from xml file
                _ecpbasisset.LoadPseudopotentialSet(_ecp_name);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded ECP library " << _ecp_name << flush;

                // fill auxiliary ECP basis by going through all atoms
                _ecp.ECPFill(&_ecpbasisset, _atoms);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled ECP Basis of size " << _ecp.getNumofShells() << flush;
            }

            // setup numerical integration grid
            _gridIntegration.GridSetup(_grid_name, &_dftbasisset, _atoms, &_dftbasis);
            _gridIntegration.setXCfunctional(_xc_functional_name);

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup numerical integration grid " << _grid_name << " for vxc functional "
                    << _xc_functional_name << " with " << _gridIntegration.getGridSize() << " points" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << " divided into " << _gridIntegration.getBoxesSize() << " boxes" << flush;
            if (_use_small_grid) {
                _gridIntegration_small.GridSetup(_grid_name_small, &_dftbasisset, _atoms, &_dftbasis);
                _gridIntegration_small.setXCfunctional(_xc_functional_name);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup small numerical integration grid " << _grid_name_small << " for vxc functional "
                        << _xc_functional_name << " with " << _gridIntegration_small.getGridpoints().size() << " points" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << " divided into " << _gridIntegration_small.getBoxesSize() << " boxes" << flush;
            }

            if (_do_externalfield) {
                _gridIntegration_ext.GridSetup(_grid_name_ext, &_dftbasisset, _atoms, &_dftbasis);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup numerical integration grid " << _grid_name_ext
                        << " for external field with " << _gridIntegration_ext.getGridpoints().size() << " points" << flush;
            }

            Elements _elements;
            //set number of electrons and such


            for (unsigned i = 0; i < _atoms.size(); i++) {
                _numofelectrons += _elements.getNucCrg(_atoms[i]->type);
            }

            // if ECP
            if (_with_ecp) {
                for (unsigned i = 0; i < _atoms.size(); i++) {
                    if (_atoms[i]->type == "H" || _atoms[i]->type == "He") {
                        continue;
                    } else {
                        _numofelectrons -= _ecpbasisset.getElement(_atoms[i]->type)->getNcore();
                    }
                }
            }
            // here number of electrons is actually the total number, everywhere else in votca it is just alpha_electrons

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Total number of electrons: " << _numofelectrons << flush;

            ConfigOrbfile(_orbitals);
            }
            SetupInvariantMatrices();
            return;
        }

        ub::matrix<double> DFTENGINE::DensityMatrix_unres(const ub::matrix<double>& MOs, int numofelec) {
            if (numofelec == 0) {
                return ub::zero_matrix<double>(MOs.size1());
            }

            ub::matrix<double> _dmatGS = ub::zero_matrix<double>(MOs.size1());
#pragma omp parallel for
            for (unsigned _i = 0; _i < MOs.size1(); _i++) {
                for (unsigned _j = 0; _j < MOs.size1(); _j++) {
                    for (int _level = 0; _level < numofelec; _level++) {

                        _dmatGS(_i, _j) += MOs(_level, _i) * MOs(_level, _j);

                    }
                }
            }
            //}
            // return
            return _dmatGS;
        }

        ub::matrix<double> DFTENGINE::DensityMatrix_frac(const ub::matrix<double>& MOs, const ub::vector<double>& MOEnergies, int numofelec) {
            if (numofelec == 0) {
                return ub::zero_matrix<double>(MOs.size1());
            }

            ub::vector<double>occupation = ub::zero_vector<double>(MOEnergies.size());

            double buffer = 0.0001;
            double homo_energy = MOEnergies(numofelec - 1);
            std::vector<unsigned> degeneracies;

            for (unsigned _level = 0; _level < occupation.size(); _level++) {
                if (MOEnergies(_level)<(homo_energy - buffer)) {
                    occupation(_level) = 1.0;
                    numofelec--;
                } else if (std::abs(MOEnergies(_level) - homo_energy) < buffer) {
                    degeneracies.push_back(_level);
                } else if (MOEnergies(_level)>(homo_energy + buffer)) {
                    occupation(_level) = 0.0;
                }
            }
            double deg_occupation = double(numofelec) / double(degeneracies.size());
            for (unsigned _level = 0; _level < degeneracies.size(); _level++) {
                occupation(degeneracies[_level]) = deg_occupation;
            }
            ub::matrix<double> _dmatGS = ub::zero_matrix<double>(MOs.size1());
#pragma omp parallel for
            for (unsigned _i = 0; _i < MOs.size1(); _i++) {
                for (unsigned _j = 0; _j < MOs.size1(); _j++) {
                    for (unsigned _level = 0; _level < occupation.size(); _level++) {

                        _dmatGS(_i, _j) += occupation(_level) * MOs(_level, _i) * MOs(_level, _j);

                    }
                }
            }
            return _dmatGS;
        }

        void DFTENGINE::NuclearRepulsion() {
            Elements element;
            E_nucnuc = 0.0;

            std::vector<double> charge;
            for (unsigned i = 0; i < _atoms.size(); i++) {
                string name = _atoms[i]->type;
                //cout << " Using atom " << name << "\n" << endl;
                double Q = element.getNucCrg(name);
                bool HorHe = (name == "H" || name == "He");
                if (_with_ecp && !HorHe) {
                    Q -= _ecpbasisset.getElement(name)->getNcore();
                }
                charge.push_back(Q);
            }


            for (unsigned i = 0; i < _atoms.size(); i++) {
                const tools::vec& r1 = _atoms[i]->getPos() * tools::conv::ang2bohr;
                double charge1 = charge[i];
                for (unsigned j = 0; j < i; j++) {
                    const tools::vec& r2 = _atoms[j]->getPos() * tools::conv::ang2bohr;
                    double charge2 = charge[j];
                    E_nucnuc += charge1 * charge2 / (abs(r1 - r2));
                }
            }
            return;
        }

        double DFTENGINE::ExternalRepulsion(ctp::Topology* top) {
            Elements element;

            if (_externalsites.size() == 0) {
                return 0;
            }

            QMMInterface qmminter;
            ctp::PolarSeg nuclei = qmminter.Convert(_atoms);

            ctp::PolarSeg::iterator pes;
            for (ctp::APolarSite* nucleus:nuclei) {
                nucleus->setIsoP(0.0);
                string name = nucleus->getName();
                double Q = element.getNucCrg(name);
                bool HorHe = (name == "H" || name == "He");
                if (_with_ecp && !HorHe) {
                    Q -= _ecpbasisset.getElement(name)->getNcore();
                }
                nucleus->setQ00(Q, 0);
            }
            ctp::XInteractor actor;
            actor.ResetEnergy();
            nuclei.CalcPos();
            double E_ext = 0.0;
            for (ctp::PolarSeg* seg:_externalsites) {
                seg->CalcPos();
                tools::vec s = tools::vec(0.0);
                if (top == NULL) {
                    s = nuclei.getPos() -seg->getPos();
                } else {
                    s = top->PbShortestConnect(nuclei.getPos(),seg->getPos()) + nuclei.getPos() - seg->getPos();
                }

                for (auto nucleus:nuclei) {
                    for (auto site:(*seg)) {
                        actor.BiasIndu(*nucleus, *site, s);
                        nucleus->Depolarize();
                        E_ext += actor.E_f(*nucleus, *site);


                    }
                }
            }
            return E_ext * tools::conv::int2eV * tools::conv::ev2hrt;
        }

        double DFTENGINE::ExternalGridRepulsion(std::vector<double> externalpotential_nuc) {
            Elements element;
            double E_ext = 0.0;

            if (!_do_externalfield) {
                return 0;
            }

            for (unsigned i = 0; i < _atoms.size(); i++) {
                string name = _atoms[i]->type;
                double Q = element.getNucCrg(name);
                bool HorHe = (name == "H" || name == "He");
                if (_with_ecp && !HorHe) {
                    Q -= _ecpbasisset.getElement(name)->getNcore();
                }
                E_ext += Q * externalpotential_nuc[i];
            }

            return E_ext;
        }

        string DFTENGINE::Choosesmallgrid(string largegrid) {
            string smallgrid;

            if (largegrid == "xfine") {
                smallgrid = "fine";
            } else if (largegrid == "fine") {
                smallgrid = "medium";
            } else if (largegrid == "medium") {
                _use_small_grid = false;
                smallgrid = "medium";
            } else if (largegrid == "coarse") {
                _use_small_grid = false;
                smallgrid = "coarse";
            } else if (largegrid == "xcoarse") {
                _use_small_grid = false;
                smallgrid = "xcoarse";
            } else {
                throw runtime_error("Grid name for Vxc integration not known.");
            }

            return smallgrid;
        }


        //average atom densities matrices, for SP and other combined shells average each subshell separately. Does not really work yet!!

        ub::matrix<double> DFTENGINE::AverageShells(const ub::matrix<double>& dmat, AOBasis& dftbasis) {
            ub::matrix<double> avdmat = ub::zero_matrix<double>(dmat.size1());
            AOBasis::AOShellIterator it;
            int start = 0.0;
            std::vector<int> starts;
            std::vector<int> ends;
            for (it = dftbasis.firstShell(); it < dftbasis.lastShell(); ++it) {
                const AOShell* shell = dftbasis.getShell(it);
                int end = shell->getNumFunc() + start;

                if (shell->getLmax() != shell->getLmin()) {
                    std::vector<int> temp = NumFuncSubShell(shell->getType());
                    int numfunc = start;
                    for (unsigned i = 0; i < temp.size(); i++) {

                        starts.push_back(numfunc);
                        numfunc += temp[i];
                        ends.push_back(numfunc);
                    }
                } else {
                    starts.push_back(start);
                    ends.push_back(end);
                }
                start = end;
            }
            for (unsigned k = 0; k < starts.size(); k++) {
                int s1 = starts[k];
                int e1 = ends[k];
                int len1 = e1 - s1;
                for (unsigned l = 0; l < starts.size(); l++) {
                    int s2 = starts[l];
                    int e2 = ends[l];
                    int len2 = e2 - s2;
                    double diag = 0.0;
                    double offdiag = 0.0;
                    for (int i = 0; i < len1; ++i) {
                        for (int j = 0; j < len2; ++j) {
                            if (i == j) {
                                diag += dmat(s1 + i, s2 + j);
                            } else {
                                offdiag += dmat(s1 + i, s2 + j);
                            }
                        }
                    }
                    if (len1 == len2) {
                        diag = diag / double(len1);
                        offdiag = offdiag / double(len1 * (len1 - 1));
                    } else {
                        double avg = (diag + offdiag) / double(len1 * len2);
                        diag = avg;
                        offdiag = avg;
                    }
                    for (int i = 0; i < len1; ++i) {
                        for (int j = 0; j < len2; ++j) {
                            if (i == j) {
                                avdmat(s1 + i, s2 + j) = diag;
                            } else {
                                avdmat(s1 + i, s2 + j) = offdiag;
                            }
                        }
                    }
                }
            }

            return avdmat;
        }
    }
}
