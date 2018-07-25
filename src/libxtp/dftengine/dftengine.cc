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

// Overload of uBLAS prod function with MKL/GSL implementations

#include "votca/xtp/aobasis.h"
#include "votca/xtp/qminterface.h"
#include "votca/xtp/qmatom.h"
#include <votca/xtp/dftengine.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/xtp/aomatrix.h>

#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <votca/tools/constants.h>

#include <votca/tools/elements.h>

#include <votca/ctp/xinteractor.h>
#include <votca/ctp/logger.h>



using boost::format;
using namespace boost::filesystem;
using std::flush;
using namespace votca::tools;

namespace votca {
  namespace xtp {

    // +++++++++++++++++++++++++++++ //
    // DFTENGINE MEMBER FUNCTIONS        //
    // +++++++++++++++++++++++++++++ //

    void DFTENGINE::CleanUp() {

    }

    void DFTENGINE::Initialize(Property& options) {

      string key = "package";

      // get OpenMP thread number

      _openmp_threads = options.ifExistsReturnElseReturnDefault<int>(key + ".openmp", 0);

      // basis sets
      _dftbasis_name = options.ifExistsReturnElseThrowRuntimeError<string>(key + ".dftbasis");

      if (options.exists(key + ".auxbasis")) {
        _auxbasis_name = options.get(key + ".auxbasis").as<string>();
        _with_RI = true;
      } else {
        _with_RI = false;
      }

      if (!_with_RI) {
        std::vector<std::string> choices = {"direct", "cache"};
        _four_center_method = options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(key + ".four_center_method", choices);

        _with_screening = options.ifExistsReturnElseReturnDefault<bool>(key + ".with_screening", true);
        _screening_eps = options.ifExistsReturnElseReturnDefault<double>(key + ".screening_eps", 1e-9);
      }

      if (options.exists(key + ".ecp")) {
        _ecp_name = options.get(key + ".ecp").as<string>();
        _with_ecp = true;
      } else {
        _with_ecp = false;
      }
      _with_guess = options.ifExistsReturnElseReturnDefault<bool>(key + ".read_guess", false);
      _initial_guess = options.ifExistsReturnElseReturnDefault<string>(key + ".initial_guess", "atom");


      // numerical integrations
      _grid_name = options.ifExistsReturnElseReturnDefault<string>(key + ".integration_grid", "medium");
      _use_small_grid = options.ifExistsReturnElseReturnDefault<bool>(key + ".integration_grid_small", true);
      _grid_name_small = Choosesmallgrid(_grid_name);

      // exchange and correlation as in libXC

      _xc_functional_name = options.ifExistsReturnElseThrowRuntimeError<string>(key + ".xc_functional");

      _numofelectrons = 0;

      if (options.exists(key + ".convergence")) {

        _Econverged = options.ifExistsReturnElseReturnDefault<double>(key + ".convergence.energy", 1e-7);
        _error_converged = options.ifExistsReturnElseReturnDefault<double>(key + ".convergence.error", 1e-7);
        _max_iter = options.ifExistsReturnElseReturnDefault<int>(key + ".convergence.max_iterations", 100);



        if (options.exists(key + ".convergence.method")) {
          string method = options.get(key + ".convergence.method").as<string>();
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

        _mixingparameter = options.ifExistsReturnElseReturnDefault<double>(key + ".convergence.mixing", -10.0);
        _levelshift = options.ifExistsReturnElseReturnDefault<double>(key + ".convergence.levelshift", 0.0);
        _levelshiftend = options.ifExistsReturnElseReturnDefault<double>(key + ".convergence.levelshift_end", 0.8);
        _maxout = options.ifExistsReturnElseReturnDefault<bool>(key + ".convergence.DIIS_maxout", false);
        _histlength = options.ifExistsReturnElseReturnDefault<int>(key + ".convergence.DIIS_length", 10);
        _diis_start = options.ifExistsReturnElseReturnDefault<double>(key + ".convergence.DIIS_start", 0.01);
        _adiis_start = options.ifExistsReturnElseReturnDefault<double>(key + ".convergence.ADIIS_start", 2);

      } else {
        _Econverged = 1e-7;
        _error_converged = 1e-7;
        _maxout = false;
        _diis_start = 0.01;
        _adiis_start = 2;
        _histlength = 10;
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




    void DFTENGINE::PrintMOs(const Eigen::VectorXd& MOEnergies){
      CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t Orbital energies: " << flush;
      CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t index occupation energy(Hartree) " << flush;
      for (int i = 0; i<MOEnergies.size(); i++) {
        int occupancy=0;
        if (i < _numofelectrons / 2 ) {
          occupancy=2;
        }
          CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << i << "\t"<<occupancy<<"\t"<< std::setprecision(12) << MOEnergies(i) << flush;
      }
      return;
    }

    bool DFTENGINE::Evaluate(Orbitals& orbitals) {

#ifdef _OPENMP
      omp_set_num_threads(_openmp_threads);
#endif

      Eigen::VectorXd& MOEnergies = orbitals.MOEnergies();
      Eigen::MatrixXd& MOCoeff = orbitals.MOCoefficients();
      if (MOEnergies.size() != _dftbasis.AOBasisSize()) {
        MOEnergies.resize(_dftbasis.AOBasisSize());
      }
      if (MOCoeff.rows() != _dftbasis.AOBasisSize() || MOCoeff.cols() != _dftbasis.AOBasisSize()) {
        MOCoeff.conservativeResize(_dftbasis.AOBasisSize(), _dftbasis.AOBasisSize());
      }

      Eigen::MatrixXd H0 = _dftAOkinetic.Matrix() + _dftAOESP.getNuclearpotential();
      NuclearRepulsion();

      if(_with_ecp){
          H0+=_dftAOECP.Matrix();
      }

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
        double externalgrid_nucint = ExternalGridRepulsion(_externalgrid_nuc);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Nuclei external potential interaction " << externalgrid_nucint << " Hartree" << flush;
        H0 += _gridIntegration_ext.IntegrateExternalPotential(_externalgrid);
        E_nucnuc += externalgrid_nucint;
      }

      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Nuclear Repulsion Energy is " << E_nucnuc << flush;


      // if we have a guess we do not need this.
      if (_with_guess) {
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Reading guess from orbitals object/file" << flush;
        _dftAOdmat = orbitals.DensityMatrixGroundState();
      } else if (guess_set) {
        ConfigOrbfile(orbitals);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using starting guess from last QM/MM iteration" << flush;
        _dftAOdmat = last_dmat;
      } else {
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup Initial Guess using: " << _initial_guess << flush;
        if (_initial_guess == "independent") {
          conv_accelerator.SolveFockmatrix(MOEnergies, MOCoeff, H0);
          _dftAOdmat = orbitals.DensityMatrixGroundState();

        } else if (_initial_guess == "atom") {

          _dftAOdmat = AtomicGuess(orbitals);
          CalculateERIs(_dftbasis, _dftAOdmat);

          if (_use_small_grid) {
            orbitals.AOVxc() = _gridIntegration_small.IntegrateVXC(_dftAOdmat);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled approximate DFT Vxc matrix " << flush;
          } else {
            orbitals.AOVxc() = _gridIntegration.IntegrateVXC(_dftAOdmat);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Vxc matrix " << flush;
          }
          Eigen::MatrixXd H = H0 + _ERIs.getERIs() + orbitals.AOVxc();
          if(_ScaHFX>0){
              if (_with_RI) {
             _ERIs.CalculateEXX(_dftAOdmat);
           } else {
             _ERIs.CalculateEXX_4c_small_molecule(_dftAOdmat);
           }
            H-=0.5*_ScaHFX*_ERIs.getEXX();
          }
          
          conv_accelerator.SolveFockmatrix(MOEnergies, MOCoeff, H);
          _dftAOdmat = orbitals.DensityMatrixGroundState();

          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Full atomic density Matrix gives N=" << std::setprecision(9) << _dftAOdmat.cwiseProduct(_dftAOoverlap.Matrix()).sum() << " electrons." << flush;
        } else {
          throw runtime_error("Initial guess method not known/implemented");
        }

      }

      orbitals.setQMpackage("xtp");

      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " STARTING SCF cycle" << flush;
      CTP_LOG(ctp::logDEBUG, *_pLog) << " --------------------------------------------------------------------------" << flush;

      double energyold = std::numeric_limits<double>::max();

      for (int _this_iter = 0; _this_iter < _max_iter; _this_iter++) {
        CTP_LOG(ctp::logDEBUG, *_pLog) << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Iteration " << _this_iter + 1 << " of " << _max_iter << flush;

        CalculateERIs(_dftbasis, _dftAOdmat);

        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Electron repulsion matrix of dimension: " << _ERIs.getSize1() << " x " << _ERIs.getSize2() << flush;
        double vxcenergy = 0.0;
        if (_use_small_grid && conv_accelerator.getDIIsError() > 1e-3) {
          orbitals.AOVxc() = _gridIntegration_small.IntegrateVXC(_dftAOdmat);
          vxcenergy = _gridIntegration_small.getTotEcontribution();
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled approximate DFT Vxc matrix " << flush;
        } else {
          orbitals.AOVxc() = _gridIntegration.IntegrateVXC(_dftAOdmat);
          vxcenergy = _gridIntegration.getTotEcontribution();
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Vxc matrix " << flush;
        }
        
        Eigen::MatrixXd H = H0 + _ERIs.getERIs() + orbitals.AOVxc();
        if(_ScaHFX>0){
              if (_with_RI) {
                if(conv_accelerator.getUseMixing()){
                  _ERIs.CalculateEXX(_dftAOdmat); 
                }else{
                  Eigen::Block<Eigen::MatrixXd> occblock=MOCoeff.block(0,0,MOEnergies.rows(), _numofelectrons / 2);
                  _ERIs.CalculateEXX(occblock,_dftAOdmat); 
                }
           } else {
             _ERIs.CalculateEXX_4c_small_molecule(_dftAOdmat);
           } 
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Electron exchange matrix of dimension: " << _ERIs.getSize1() << " x " << _ERIs.getSize2() << flush;
            H-=0.5*_ScaHFX*_ERIs.getEXX();
          }
      
        double Eone = _dftAOdmat.cwiseProduct(H0).sum();
        double Etwo = 0.5 * _ERIs.getERIsenergy() + vxcenergy;
        if(_ScaHFX>0){
          Etwo-=_ScaHFX/4*_ERIs.getEXXsenergy();
        }
        double totenergy = Eone + E_nucnuc + Etwo;
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Single particle energy " << std::setprecision(12) << Eone << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Two particle energy " << std::setprecision(12) << Etwo << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << std::setprecision(12) << " Local Exc contribution " << vxcenergy << flush;
        if(_ScaHFX>0){
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << std::setprecision(12) << " Non local Ex contribution " << -_ScaHFX/4*_ERIs.getEXXsenergy() << flush;
        }
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Total Energy " << std::setprecision(12) << totenergy << flush;

        _dftAOdmat = conv_accelerator.Iterate(_dftAOdmat, H, MOEnergies, MOCoeff, totenergy);

        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " DIIs error " << conv_accelerator.getDIIsError() << std::flush;
        double deltaE=totenergy - energyold;
        if(_this_iter == 0){
          deltaE=0;
        }
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Delta Etot " <<deltaE<< std::flush;
        if (tools::globals::verbose) {
          PrintMOs(MOEnergies);
        }

        CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\tGAP " << MOEnergies(_numofelectrons / 2) - MOEnergies(_numofelectrons / 2 - 1) << flush;

        if (std::abs(totenergy - energyold) < _Econverged && conv_accelerator.getDIIsError() < _error_converged) {
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Total Energy has converged to " << std::setprecision(9) << std::abs(totenergy - energyold) << "[Ha] after " << _this_iter + 1 <<
                  " iterations. DIIS error is converged up to " << _error_converged << "[Ha]" << flush;
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Final Single Point Energy " << std::setprecision(12) << totenergy << " Ha" << flush;
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " MO Energies  [Ha]" << flush;
          
          PrintMOs(MOEnergies);

          last_dmat = _dftAOdmat;
          guess_set = true;
          // orbitals saves total energies in [eV]
          orbitals.setQMEnergy(totenergy * tools::conv::hrt2ev);
          break;
        } else {
          energyold = totenergy;
        }
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Density Matrix gives N=" << std::setprecision(9) << _dftAOdmat.cwiseProduct(_dftAOoverlap.Matrix()).sum() << " electrons." << flush;

      }
      return true;
    }




    // SETUP INVARIANT AOMatrices

    void DFTENGINE::SetupInvariantMatrices() {
      
        _dftAOoverlap.Fill(_dftbasis);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftAOoverlap.Dimension() << flush;

      _dftAOkinetic.Fill(_dftbasis);
      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Kinetic energy matrix of dimension: " << _dftAOkinetic.Dimension() << flush;

      _dftAOESP.Fillnucpotential(_dftbasis, _atoms);
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
            CTP_LOG(ctp::logDEBUG, *_pLog) << flush;
          }
        }
      }

      if (_with_ecp) {
        _dftAOECP.setECP(&_ecp);
        _dftAOECP.Fill(_dftbasis);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT ECP matrix of dimension: " << _dftAOECP.Dimension() << flush;
      }

      conv_accelerator.Configure(ConvergenceAcc::KSmode::closed, _usediis,
              true, _histlength, _maxout, _adiis_start, _diis_start, _levelshift, _levelshiftend, _numofelectrons / 2, _mixingparameter);
      conv_accelerator.setLogger(_pLog);
      conv_accelerator.setOverlap(&_dftAOoverlap.Matrix(), 1e-8);

      if (_with_RI) {

        AOCoulomb _auxAOcoulomb;
        _auxAOcoulomb.Fill(_auxbasis);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled AUX Coulomb matrix of dimension: " << _auxAOcoulomb.Dimension() << flush;

        Eigen::MatrixXd Inverse=_auxAOcoulomb.Pseudo_InvSqrt(1e-8);

        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Inverted AUX Coulomb matrix, removed " << _auxAOcoulomb.Removedfunctions() << " functions from aux basis" << flush;

        // prepare invariant part of electron repulsion integrals
        _ERIs.Initialize(_dftbasis, _auxbasis, Inverse);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup invariant parts of Electron Repulsion integrals " << flush;
      } else {

        if (_four_center_method=="cache") {

          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculating 4c integrals. " << flush;
          _ERIs.Initialize_4c_small_molecule(_dftbasis);
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculated 4c integrals. " << flush;
        }

        if (_with_screening && _four_center_method=="direct") {

          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculating 4c diagonals. " << flush;
          _ERIs.Initialize_4c_screening(_dftbasis, _screening_eps);
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculated 4c diagonals. " << flush;
        }
      }

      return;
    }

    Eigen::MatrixXd DFTENGINE::AtomicGuess(Orbitals& orbitals) {
      Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(_dftbasis.AOBasisSize(), _dftbasis.AOBasisSize());

      std::vector<QMAtom*> uniqueelements;
      std::vector<QMAtom*>::const_iterator at;
      std::vector<QMAtom*>::iterator st;
      std::vector< Eigen::MatrixXd > uniqueatom_guesses;
      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Scanning molecule of size " << _atoms.size() << " for unique elements" << flush;
      for (QMAtom* atom:_atoms) {
        bool exists = false;
        if (uniqueelements.size() == 0) {
          exists = false;
        } else {
          for (QMAtom* uniqueatom: uniqueelements) {
            if (atom->getType() == uniqueatom->getType()) {
              exists = true;
              break;
            }
          }
        }
        if (!exists) {
          uniqueelements.push_back(atom);
        }
      }
      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " " << uniqueelements.size() << " unique elements found" << flush;
      tools::Elements _elements;
      for (QMAtom* uniqueatom: uniqueelements) {

        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculating atom density for " << uniqueatom->getType() << flush;
        bool with_ecp = _with_ecp;
        if (uniqueatom->getType() == "H" || uniqueatom->getType() == "He") {
          with_ecp = false;
        }
        std::vector<QMAtom*> atom;
        atom.push_back(uniqueatom);

        AOBasis dftbasis;

        NumericalIntegration gridIntegration;
        dftbasis.AOBasisFill(_dftbasisset, atom);
        AOBasis ecp;
        if (with_ecp) {
          ecp.ECPFill(_ecpbasisset, atom);
        }
        gridIntegration.GridSetup(_grid_name, atom, &dftbasis);
        gridIntegration.setXCfunctional(_xc_functional_name);
        int numofelectrons = (*st)->getNuccharge();
        int alpha_e = 0;
        int beta_e = 0;

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
        ERIs ERIs_atom;

        // DFT AOOverlap matrix

        dftAOoverlap.Fill(dftbasis);
        dftAOkinetic.Fill(dftbasis);


        dftAOESP.Fillnucpotential(dftbasis, atom);
        ERIs_atom.Initialize_4c_small_molecule(dftbasis);

        Eigen::VectorXd MOEnergies_alpha;
        Eigen::MatrixXd MOCoeff_alpha;
        Eigen::VectorXd MOEnergies_beta;
        Eigen::MatrixXd MOCoeff_beta;
        ConvergenceAcc Convergence_alpha;

        ConvergenceAcc Convergence_beta;
        double adiisstart = 0;
        double diisstart = 0;

        Convergence_alpha.Configure(ConvergenceAcc::KSmode::open, true, false, 20, 0, adiisstart, diisstart, 0.1, diisstart, alpha_e, -1);
        Convergence_alpha.setLogger(_pLog);
        Convergence_alpha.setOverlap(&dftAOoverlap.Matrix(), 1e-8);
        Convergence_beta.Configure(ConvergenceAcc::KSmode::open, true, false, 20, 0, adiisstart, diisstart, 0.1, diisstart, beta_e, -1);
        Convergence_beta.setLogger(_pLog);
        Convergence_beta.setOverlap(&dftAOoverlap.Matrix(), 1e-8);
        /**** Construct initial density  ****/

        Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.getNuclearpotential();
        if (with_ecp) {
          dftAOECP.setECP(&ecp);
          dftAOECP.Fill(dftbasis);
          H0 += dftAOECP.Matrix();
        }
        Convergence_alpha.SolveFockmatrix(MOEnergies_alpha, MOCoeff_alpha, H0);

        Eigen::MatrixXd dftAOdmat_alpha = Convergence_alpha.DensityMatrix(MOCoeff_alpha, MOEnergies_alpha);
        if ((*st)->getType() == "H") {
          uniqueatom_guesses.push_back(dftAOdmat_alpha);
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Atomic density Matrix for " << (*st)->getType() <<
                  " gives N=" << std::setprecision(9) << dftAOdmat_alpha.cwiseProduct(dftAOoverlap.Matrix()).sum() << " electrons." << flush;
          continue;
        }
        Convergence_beta.SolveFockmatrix(MOEnergies_beta, MOCoeff_beta, H0);
        Eigen::MatrixXd dftAOdmat_beta = Convergence_beta.DensityMatrix(MOCoeff_beta, MOEnergies_beta);

        
        double energyold = 0;
        int maxiter = 50;
        for (int this_iter = 0; this_iter < maxiter; this_iter++) {
          ERIs_atom.CalculateERIs_4c_small_molecule(dftAOdmat_alpha + dftAOdmat_beta);
          double E_two_alpha = ERIs_atom.getERIs().cwiseProduct(dftAOdmat_alpha).sum();
          double E_two_beta = ERIs_atom.getERIs().cwiseProduct(dftAOdmat_beta).sum();
          Eigen::MatrixXd H_alpha = H0 + ERIs_atom.getERIs();
          Eigen::MatrixXd H_beta = H0 + ERIs_atom.getERIs();
          
            Eigen::MatrixXd AOVxc_alpha = gridIntegration.IntegrateVXC(dftAOdmat_alpha);
            double E_vxc_alpha = gridIntegration.getTotEcontribution();

            Eigen::MatrixXd AOVxc_beta = gridIntegration.IntegrateVXC(dftAOdmat_beta);
            double E_vxc_beta = gridIntegration.getTotEcontribution();
            H_alpha += AOVxc_alpha;
            H_beta += AOVxc_beta;
            E_two_alpha += E_vxc_alpha;
            E_two_beta += E_vxc_beta;
          
            if(_ScaHFX>0){
              ERIs_atom.CalculateEXX_4c_small_molecule(dftAOdmat_alpha);
              double E_exx_alpha =-0.5*_ScaHFX*ERIs_atom.getEXX().cwiseProduct(dftAOdmat_alpha).sum();
              H_alpha -=_ScaHFX* ERIs_atom.getEXX();
              E_two_alpha += E_exx_alpha;
              ERIs_atom.CalculateEXX_4c_small_molecule(dftAOdmat_beta);
              double E_exx_beta =-0.5*_ScaHFX*ERIs_atom.getEXX().cwiseProduct(dftAOdmat_beta).sum();
              H_beta -=_ScaHFX*ERIs_atom.getEXX();
              E_two_beta += E_exx_beta;
            }

          double E_one_alpha = dftAOdmat_alpha.cwiseProduct(H0).sum();
          double E_one_beta = dftAOdmat_beta.cwiseProduct(H0).sum();
          double E_alpha = E_one_alpha + E_two_alpha;
          double E_beta = E_one_beta + E_two_beta;
          double totenergy = E_alpha + E_beta;
          //evolve alpha
          dftAOdmat_alpha = Convergence_alpha.Iterate(dftAOdmat_alpha, H_alpha, MOEnergies_alpha, MOCoeff_alpha, E_alpha);
          //evolve beta
          dftAOdmat_beta = Convergence_beta.Iterate(dftAOdmat_beta, H_beta, MOEnergies_beta, MOCoeff_beta, E_beta);

          if (tools::globals::verbose) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Iter " << this_iter << " of " << maxiter
                    << " Etot " << totenergy << " diise_a " << Convergence_alpha.getDIIsError() << " diise_b " << Convergence_beta.getDIIsError()
                    << " a_gap " << MOEnergies_alpha(alpha_e) - MOEnergies_alpha(alpha_e - 1) << " b_gap " << MOEnergies_beta(beta_e) - MOEnergies_beta(beta_e - 1) << flush;
          }
          bool converged = (std::abs(totenergy - energyold) < _Econverged &&
                  Convergence_alpha.getDIIsError() < _error_converged &&
                  Convergence_beta.getDIIsError() < _error_converged);
          if (converged || this_iter == maxiter - 1) {

            if (converged) {
              CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Converged after " << this_iter + 1 << " iterations" << flush;
            } else {
              CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Not converged after " << this_iter + 1 <<
                      " iterations. Using unconverged density. DIIsError_alpha=" << Convergence_alpha.getDIIsError()
                      << " DIIsError_beta=" << Convergence_beta.getDIIsError() << flush;
            }


            Eigen::MatrixXd avdmat = AverageShells(dftAOdmat_alpha + dftAOdmat_beta, dftbasis);
            uniqueatom_guesses.push_back(avdmat);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Atomic density Matrix for " << (*st)->getType() << " gives N="
                    << std::setprecision(9) << avdmat.cwiseProduct(dftAOoverlap.Matrix()).sum() << " electrons." << flush;
            break;

          } else {
            energyold = totenergy;
          }
        }
      }
      unsigned start = 0;


      for (at = _atoms.begin(); at < _atoms.end(); ++at) {
        unsigned index = 0;
        for (unsigned i = 0; i < uniqueelements.size(); i++) {
          if ((*at)->getType() == uniqueelements[i]->getType()) {
            index = i;
            break;
          }
        }
        double size = _dftbasis.getFuncperAtom((*at)->getAtomID());
        guess.block(start, start, size, size) = uniqueatom_guesses[index];

        start += size;
      }

      return guess;
    }

    void DFTENGINE::ConfigOrbfile(Orbitals& orbitals) {
      if (_with_guess) {

        if (orbitals.hasDFTbasis()) {
          if (orbitals.getDFTbasis() != _dftbasis_name) {
            throw runtime_error((boost::format("Basisset Name in guess orb file and in dftengine option file differ %1% vs %2%") % orbitals.getDFTbasis() % _dftbasis_name).str());
          }
        } else {
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " WARNING: Orbital file has no basisset information,using it as a guess might work or not for calculation with " << _dftbasis_name << flush;
        }
      }
      orbitals.setDFTbasis(_dftbasis_name);
      orbitals.setBasisSetSize(_dftbasis.AOBasisSize());
      orbitals.setScaHFX(_ScaHFX);
      if (_with_ecp) {
        orbitals.setECP(_ecp_name);
      }
      if (_with_RI) {
        orbitals.setAuxbasis(_auxbasis_name);
      }

      if (_with_guess) {
        if (orbitals.hasECP() || _with_ecp) {
          if (orbitals.getECP() != _ecp_name) {
            throw runtime_error((boost::format("ECPs in orb file: %1% and options %2% differ") % orbitals.getECP() % _ecp_name).str());
          }
        }
        if (orbitals.getNumberOfElectrons() != _numofelectrons / 2) {
          throw runtime_error((boost::format("Number of electron in guess orb file: %1% and in dftengine: %2% differ.") % orbitals.getNumberOfElectrons() % (_numofelectrons / 2)).str());
        }
        if (orbitals.getNumberOfLevels() != _dftbasis.AOBasisSize()) {
          throw runtime_error((boost::format("Number of levels in guess orb file: %1% and in dftengine: %2% differ.") % orbitals.getNumberOfLevels() % _dftbasis.AOBasisSize()).str());
        }
      } else {
        orbitals.setNumberOfElectrons(_numofelectrons / 2);
        orbitals.setNumberOfLevels(_numofelectrons / 2, _dftbasis.AOBasisSize() - _numofelectrons / 2);
      }
      return;
    }


    // PREPARATION

    void DFTENGINE::Prepare(Orbitals& orbitals) {
#ifdef _OPENMP

      omp_set_num_threads(_openmp_threads);
      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using " << omp_get_max_threads() << " threads" << flush;

#endif
      if(tools::globals::VOTCA_MKL){
     CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()
                                 << " Using MKL overload for Eigen "<< flush;
  }else{
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()
                                 << " Using native Eigen implementation, no BLAS overload " << flush;
  }

      if (_atoms.size() == 0) {
        _atoms = orbitals.QMAtoms();
      }

      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Molecule Coordinates [A] " << flush;
      for (unsigned i = 0; i < _atoms.size(); i++) {
        CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << _atoms[i]->getType() << " " << _atoms[i]->getPos() * conv::bohr2ang << " " << flush;
      }

      // load and fill DFT basis set
      _dftbasisset.LoadBasisSet(_dftbasis_name);

      _dftbasis.AOBasisFill(_dftbasisset, _atoms);
      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << flush;

      if (_with_RI) {
        // load and fill AUX basis set
        _auxbasisset.LoadBasisSet(_auxbasis_name);
        //_orbitals.setDFTbasis( _dftbasis_name );
        _auxbasis.AOBasisFill(_auxbasisset, _atoms);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded AUX Basis Set " << _auxbasis_name << flush;
      }
      if (_with_ecp) {
        // load ECP (element-wise information) from xml file
        _ecpbasisset.LoadPseudopotentialSet(_ecp_name);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded ECP library " << _ecp_name << flush;

        // fill auxiliary ECP basis by going through all atoms
        _ecp.ECPFill(_ecpbasisset, _atoms);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled ECP Basis of size " << _ecp.getNumofShells() << flush;
      }

      // setup numerical integration grid
      _gridIntegration.GridSetup(_grid_name, _atoms, &_dftbasis);
      _gridIntegration.setXCfunctional(_xc_functional_name);
      
      _ScaHFX=_gridIntegration.getExactExchange(_xc_functional_name);
      if(_ScaHFX>0){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using hybrid functional with alpha="<<_ScaHFX<< flush;
      }

      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup numerical integration grid " << _grid_name << " for vxc functional "
              << _xc_functional_name << " with " << _gridIntegration.getGridSize() << " points" << flush;
      CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << " divided into " << _gridIntegration.getBoxesSize() << " boxes" << flush;
      if (_use_small_grid) {
        _gridIntegration_small.GridSetup(_grid_name_small, _atoms, &_dftbasis);
        _gridIntegration_small.setXCfunctional(_xc_functional_name);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup small numerical integration grid " << _grid_name_small << " for vxc functional "
                << _xc_functional_name << " with " << _gridIntegration_small.getGridpoints().size() << " points" << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << "\t\t " << " divided into " << _gridIntegration_small.getBoxesSize() << " boxes" << flush;
      }

      if (_do_externalfield) {
        _gridIntegration_ext.GridSetup(_grid_name_ext, _atoms, &_dftbasis);
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup numerical integration grid " << _grid_name_ext
                << " for external field with " << _gridIntegration_ext.getGridpoints().size() << " points" << flush;
      }

      for (auto& atom : _atoms) {
        _numofelectrons += atom->getNuccharge();
      }


      // here number of electrons is actually the total number, everywhere else in votca it is just alpha_electrons

      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Total number of electrons: " << _numofelectrons << flush;

      ConfigOrbfile(orbitals);
      SetupInvariantMatrices();
      return;
    }

    void DFTENGINE::NuclearRepulsion() {
      E_nucnuc = 0.0;

      for (unsigned i = 0; i < _atoms.size(); i++) {
        const tools::vec& r1 = _atoms[i]->getPos();
        double charge1 = _atoms[i]->getNuccharge();
        for (unsigned j = 0; j < i; j++) {
          const tools::vec& r2 = _atoms[j]->getPos();
          double charge2 = _atoms[j]->getNuccharge();
          E_nucnuc += charge1 * charge2 / (abs(r1 - r2));
        }
      }
      return;
    }

    double DFTENGINE::ExternalRepulsion(ctp::Topology* top) {


      if (_externalsites.size() == 0) {
        return 0;
      }

      QMInterface qmminter;
      ctp::PolarSeg nuclei = qmminter.Convert(_atoms);

      for (unsigned i = 0; i < nuclei.size(); ++i) {
        ctp::APolarSite* nucleus = nuclei[i];
        nucleus->setIsoP(0.0);
        double Q = _atoms[i]->getNuccharge();
        nucleus->setQ00(Q, 0);
      }
      ctp::XInteractor actor;
      actor.ResetEnergy();
      nuclei.CalcPos();
      double E_ext = 0.0;
      for (ctp::PolarSeg* seg : _externalsites) {
        seg->CalcPos();
        tools::vec s = tools::vec(0.0);
        if (top == NULL) {
          s = nuclei.getPos() - seg->getPos();
        } else {
          s = top->PbShortestConnect(nuclei.getPos(), seg->getPos()) + nuclei.getPos() - seg->getPos();
        }

        for (auto nucleus : nuclei) {
          for (auto site : (*seg)) {
            actor.BiasIndu(*nucleus, *site, s);
            nucleus->Depolarize();
            E_ext += actor.E_f(*nucleus, *site);


          }
        }
      }
      return E_ext * tools::conv::int2eV * tools::conv::ev2hrt;
    }

    double DFTENGINE::ExternalGridRepulsion(std::vector<double> externalpotential_nuc) {
      double E_ext = 0.0;
      if (!_do_externalfield) {
        return 0;
      }
      for (unsigned i = 0; i < _atoms.size(); i++) {
        double Q = _atoms[i]->getNuccharge();
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


    //average atom densities matrices, for SP and other combined shells average each subshell separately.

    Eigen::MatrixXd DFTENGINE::AverageShells(const Eigen::MatrixXd& dmat, AOBasis& dftbasis) {
      Eigen::MatrixXd avdmat = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
      int start = 0.0;
      std::vector<int> starts;
      std::vector<int> ends;
      for ( const AOShell* shell:dftbasis) {
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

    void DFTENGINE::CalculateERIs(const AOBasis& dftbasis, const Eigen::MatrixXd& DMAT) {

      if (_with_RI)
        _ERIs.CalculateERIs(_dftAOdmat);
      else if (_four_center_method.compare("cache") == 0)
        _ERIs.CalculateERIs_4c_small_molecule(_dftAOdmat);
      else if (_four_center_method.compare("direct") == 0)
        _ERIs.CalculateERIs_4c_direct(_dftbasis, _dftAOdmat);
    }
  }
}
