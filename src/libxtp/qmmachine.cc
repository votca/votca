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


#include <votca/xtp/qmmachine.h>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/elements.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/aomatrix.h>


using boost::format;

namespace votca {
  namespace xtp {


    QMMachine::QMMachine(ctp::XJob *job, ctp::XInductor *xind, QMPackage *qmpack,
            Property *opt, string sfx, int nst, bool mav)
    : _job(job), _xind(xind), _qmpack(qmpack), _subthreads(nst),
    _isConverged(false) {

      string key = sfx + ".qmmmconvg";
      _crit_dR = opt->ifExistsReturnElseReturnDefault<double>(key + ".dR", 0.01); //nm
      _crit_dQ = opt->ifExistsReturnElseReturnDefault<double>(key + ".dQ", 0.01); //e
      _crit_dE_QM = opt->ifExistsReturnElseReturnDefault<double>(key + ".dQdE_QM", 0.001); //eV
      _crit_dE_MM = opt->ifExistsReturnElseReturnDefault<double>(key + ".dE_MM", _crit_dE_QM); //eV
      _maxIter = opt->ifExistsReturnElseReturnDefault<int>(key + ".max_iter", 32);
      _alpha = opt->ifExistsReturnElseReturnDefault<double>(key + ".mixing", 0.0);

      key = sfx;
      bool split_dpl = opt->ifExistsReturnElseReturnDefault<bool>(key + ".split_dpl", true);
      double dpl_spacing = opt->ifExistsReturnElseReturnDefault<double>(key + ".dpl_spacing", 1e-3);

      // check for archiving
      std::string _archiving_string = opt->ifExistsReturnElseReturnDefault<std::string>(key + ".archiving", "");
      if (_archiving_string.find("iterations") != std::string::npos) {
        _do_archive = true;
      } else {
        _do_archive = false;
      }

      // check for static or polarized qmmm
      key = sfx + ".tholemodel";
      _static_qmmm = true;
      qmpack->setWithPolarization(false);
      if (opt->exists(key + ".induce")) {
        bool induce = opt->get(key + ".induce").as<bool>();
        if (induce) {
          qmpack->setWithPolarization(true);
          qmpack->setDipoleSpacing(dpl_spacing);
        }
        _static_qmmm = !induce;
      }

      // GDMA options
      key = sfx + ".gdma";
      if (opt->exists(key)) {
        _do_gdma = true;
        string gdma_xml = opt->get(key).as<string> ();
        load_property_from_xml(_gdma_options, gdma_xml.c_str());
      } else {
        _do_gdma = false;
      }

      // GWBSE options
      key = sfx + ".gwbse";
      if (opt->exists(key)) {
        _do_gwbse = true;

        // PUT IN some useful type definitions for excited state
        // - singlet -> singlet exciton
        // - triplet -> triplet exciton
        // - quasiparticle -> GW quasiparticle (diagqp)
        // - kohn-sham -> DFT MO

        string _gwbse_xml = opt->get(key + ".gwbse_options").as<string> ();
        load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());
        _state = opt->get(key + ".state").as< int >();
        std::vector<std::string> choices={"singlet","triplet","quasiparticle"};
        _type =opt->ifExistsAndinListReturnElseThrowRuntimeError<string>(key + ".type",choices);
     
        key = sfx + ".gwbse.filter";
        if (opt->exists(key + ".oscillator_strength") && _type == "singlet") {
          _has_osc_filter = true;
          _osc_threshold = opt->get(key + ".oscillator_strength").as<double> ();
        }
        if (opt->exists(key + ".overlap")) {
          _has_overlap_filter = true;
        }
        if (opt->exists(key + ".localisation")) {
          _has_loc_filter = true;
          string temp = opt->get(key + ".localisation").as<string> ();
          Tokenizer tok_cleanup(temp, ", \n\t");
          std::vector <std::string> strings_vec;
          tok_cleanup.ToVector(strings_vec);
          if (strings_vec.size() != 2) {
            throw runtime_error("qmmmachine: Fragment and localisation threshold are not separated");
          }
          if (strings_vec[0] == "a" || strings_vec[0] == "A") {
            _localiseonA = true;
          } else if (strings_vec[0] == "b" || strings_vec[0] == "B") {
            _localiseonA = false;
          } else {
            throw runtime_error("qmmmachine: Fragment label not known, either A or B");
          }
          _loc_threshold = boost::lexical_cast<double>(strings_vec[1]);
        }

        if (opt->exists(key + ".charge_transfer")) {
          _has_dQ_filter = true;
          _dQ_threshold = opt->get(key + ".charge_transfer").as<double> ();
        }
        if (_has_dQ_filter && _has_loc_filter) {
          throw runtime_error("Cannot use localisation and charge_transfer filter at the same time.");
        }
        if (_has_osc_filter && _has_dQ_filter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- WARNING: filtering for optically active CT transition - might not make sense... " << flush;
          }
      } else {
        _do_gwbse = false;
      }

      return;
    }


    QMMachine::~QMMachine() {
      for (QMMIter* qit : _iters) {
        delete qit;
      }
      _iters.clear();
    }


    int QMMachine::Evaluate(ctp::XJob *job) {

      CTP_LOG(ctp::logINFO, *_log)
              << format("... dR %1$1.4f dQ %2$1.4f QM %3$1.4f MM %4$1.4f IT %5$d")
              % _crit_dR % _crit_dQ % _crit_dE_QM % _crit_dE_MM % _maxIter << flush;

      // FIGURE OUT CHARGE + MULTIPLICITY
      double dQ = 0.0;
      for (unsigned int i = 0; i < _job->getPolarTop()->QM0().size(); ++i) {
        dQ += _job->getPolarTop()->QM0()[i]->CalcTotQ();
      }
      int chrg = std::round(dQ);
      int spin = ((chrg < 0) ? -chrg : chrg) % 2 + 1;
      CTP_LOG(ctp::logINFO, *_log) << "... Q = " << chrg << ", 2S+1 = " << spin << flush;


      // PREPARE JOB DIRECTORY
      string jobFolder = "xjob_" + boost::lexical_cast<string>(_job->getId())
              + "_" + _job->getTag();
      bool created = boost::filesystem::create_directory(jobFolder);
      if (created) {
        CTP_LOG(ctp::logINFO, *_log) << "Created directory " << jobFolder << flush;
      }

      _qmpack->setCharge(chrg);
      _qmpack->setSpin(spin);

      int iterCnt = 0;
      int iterMax = _maxIter;
      for (; iterCnt < iterMax; ++iterCnt) {

        // check for polarized QM/MM convergence
        int info = Iterate(jobFolder, iterCnt);
        if (info != 0) {
          CTP_LOG(ctp::logERROR, *_log)
                  << format("Iterating job failed!") << flush;
          return 1;
        }
        if (_static_qmmm) {
          _isConverged = true;
          break;
        } else if (hasConverged()) {
          break;
        }

      }

      if (iterCnt == iterMax - 1 && !_isConverged) {
        CTP_LOG(ctp::logWARNING, *_log)
                << format("Not converged within %1$d iterations.") % iterMax;
        return 2;
      }

      return 0;
    }


    void QMMachine::RunGDMA(QMMIter* thisIter, string& runFolder){
      if (_qmpack->getPackageName() != "gaussian" || _qmpack->getExecutable() != "g03") {
        throw std::runtime_error(" Invalid QMPackage! " + _type + " Gaussian 03 only!");
      } else {
        GDMA gdma;
        gdma.Initialize(_gdma_options);
        gdma.setLog(_log);
        gdma.SetRunDir(runFolder);
        CTP_LOG(ctp::logINFO, *_log) << "Running GDMA " << flush;
        gdma.WriteInputFile();
        gdma.RunExternal();
        gdma.ParseOutputFile();
        thisIter->UpdateMPSFromGDMA(gdma.GetMultipoles(), _job->getPolarTop()->QM0());
      } 
    }


    void QMMachine::RunGWBSE(string& runFolder){
      GWBSE gwbse = GWBSE(orb_iter_input);
      // define own logger for GW-BSE that is written into a runFolder logfile
      ctp::Logger gwbse_logger(ctp::logDEBUG);
      gwbse_logger.setMultithreading(false);
      gwbse.setLogger(&gwbse_logger);
      gwbse_logger.setPreface(ctp::logINFO, (format("\nGWBSE INF ...")).str());
      gwbse_logger.setPreface(ctp::logERROR, (format("\nGWBSE ERR ...")).str());
      gwbse_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...")).str());
      gwbse_logger.setPreface(ctp::logDEBUG, (format("\nGWBSE DBG ...")).str());
      
      // Initialize with options
      gwbse.Initialize(_gwbse_options);
      // actual GW-BSE run
      gwbse.Evaluate();
      
      // write logger to log file
      ofstream ofs;
      string gwbse_logfile = runFolder + "/gwbse.log";
      ofs.open(gwbse_logfile.c_str(), ofstream::out);
      if (!ofs.is_open()) {
        throw runtime_error("Bad file handle: " + gwbse_logfile);
      }
      ofs << gwbse_logger << endl;
      ofs.close();
    }


    bool QMMachine::RunDFT(string& runFolder, std::vector<std::shared_ptr<ctp::PolarSeg> >& MultipolesBackground){
      _qmpack->setMultipoleBackground(MultipolesBackground);
      
      // setting RUNDIR for the external QMPackages, dummy for internal
      _qmpack->setRunDir(runFolder);
      
      CTP_LOG(ctp::logDEBUG, *_log) << "Writing input file " << runFolder << flush;
      _qmpack->WriteInputFile(orb_iter_input);
      
      orb_iter_input.WriteXYZ(runFolder + "/system_qm.xyz");
      
      _qmpack->Run(orb_iter_input);
      
      /* Parse information from the LOGFILE into orbitals, if
       * external QMPackage is run. Dummy for internal DFTENGINE.
       * The QMPackage's MOcoefficients are not automatically
       * parsed in DFT-only calculations and ESP fits for
       * polarized QMMM are expected in the LOGFILE.
       * IF the internal ESPFITs should be used, the MOcoefficients
       * need to be parsed too.
       */
      bool success= _qmpack->ParseLogFile(orb_iter_input);
      if(!success){
        return false;
      }
      success=_qmpack->ParseOrbitalsFile(orb_iter_input);
      _qmpack->CleanUp();
        return success;
    }

    bool QMMachine::Iterate(string jobFolder, int iterCnt) {

      // CREATE ITERATION OBJECT & SETUP RUN DIRECTORY
      QMMIter *thisIter = this->CreateNewIter();
      int iter = iterCnt;
      string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);

      bool created = boost::filesystem::create_directory(runFolder);
      if (created)
        CTP_LOG(ctp::logDEBUG, *_log) << "Created directory " << runFolder << flush;
      else
        CTP_LOG(ctp::logWARNING, *_log) << "Could not create directory " << runFolder << flush;

      // RUN CLASSICAL INDUCTION & SAVE
      _job->getPolarTop()->PrintPDB(runFolder + "/QM0_MM1_MM2.pdb");
      _xind->Evaluate(_job);

      if(!_xind->hasConverged()){
        throw std::runtime_error("Classical induction did not converge");
      }

      thisIter->setE_FM(_job->getEF00(), _job->getEF01(), _job->getEF02(),
              _job->getEF11(), _job->getEF12(), _job->getEM0(),
              _job->getEM1(), _job->getEM2(), _job->getETOT());


      /* Translate atoms in QM0() to QMAtoms in orbitals object
       * to be used in writing the QM input files for the
       * external QMPackages or directly in internal DFT engine.
       * DRAWBACK: QMAtom positions are fixed for all iterations
       *           unless the UPDATE function at the end of the
       *           iteration updates orb_iter_input (TO BE TESTED)
       */
      QMInterface qminterface;
      if (iterCnt == 0) qminterface.GenerateQMAtomsFromPolarSegs(_job->getPolarTop(), orb_iter_input);

      /* Generate list of polar segments in the MM1() and MM2()
       * region to be used in writing the background multipoles
       * in the external QMPackages or in the direct evaluation
       * in the internal DFT engine.
       */
      std::vector<std::shared_ptr<ctp::PolarSeg> > MultipolesBackground = qminterface.GenerateMultipoleList(_job->getPolarTop());

      bool success=RunDFT(runFolder, MultipolesBackground);
      if (!success) {
        return 1;
      }
       
      // GW-BSE starts here
      double energy_ex = 0.0;
      std::vector<int> state_index;

      if (_do_gwbse) {
        RunGWBSE(runFolder);

          // PROCESSING the GW-BSE result
          // - find the excited state of interest
           if (_state > 0) {
          CTP_LOG(ctp::logDEBUG, *_log) << "Excited state via GWBSE: " << flush;
          CTP_LOG(ctp::logDEBUG, *_log) << "  --- type:              " << _type << flush;
          CTP_LOG(ctp::logDEBUG, *_log) << "  --- state:             " << _state << flush;
          if (_has_overlap_filter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: overlap  " << flush;
          }
          if (_has_osc_filter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: osc.str. > " << _osc_threshold << flush;
          }
          if (_has_dQ_filter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: crg.trs. > " << _dQ_threshold << flush;
          }
          if (_has_loc_filter) {
            if (_localiseonA) {
              CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: localisation on A > " << _loc_threshold << flush;
            } else {
              CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: localisation on B > " << _loc_threshold << flush;
            }
          }
          
          

          // quasiparticle filter
          if (_has_overlap_filter) {
            if (iter == 0) {

              // One - to - One LIST in 0th iteration
              for (unsigned i = 0; i < orb_iter_input.QPdiagEnergies().size(); i++) {
                state_index.push_back(i);
              }

            } else {

              BasisSet dftbs;
              dftbs.LoadBasisSet(orb_iter_input.getDFTbasis());
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " << orb_iter_input.getDFTbasis() << flush;
              AOBasis dftbasis;
              dftbasis.AOBasisFill(dftbs, orb_iter_input.QMAtoms());
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filled DFT Basis of size " << dftbasis.AOBasisSize() << flush;

              AOOverlap dftoverlap;
              dftoverlap.Fill(dftbasis);
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << dftoverlap.Matrix().rows() << flush;

              // 'LAMBDA' matrix of the present iteration
              Eigen::MatrixXd lambda_N = orb_iter_input.CalculateQParticleAORepresentation();

              // 'LAMBDA' matrix of the previous iteration
              string runFolder_N_1 = jobFolder + "/iter_" + boost::lexical_cast<string>(iter - 1);
              string orbfile_N_1 = runFolder_N_1 + "/system.orb";
              Orbitals orbitals_N_1;
              // load the QM data from serialized orbitals object

              CTP_LOG(ctp::logDEBUG, *_log) << " Loading QM data from " << orbfile_N_1 << flush;
              orbitals_N_1.ReadFromCpt(orbfile_N_1);

              Eigen::MatrixXd lambda_N_1 = orbitals_N_1.CalculateQParticleAORepresentation();
              // calculate QP overlaps
              Eigen::MatrixXd qpoverlaps = lambda_N * dftoverlap.Matrix() * lambda_N_1.transpose();

              // filter for max absolute value (hopefully close to 1)
              for (unsigned j = 0; j < qpoverlaps.cols(); j++) {             
                int maxi = 0;
                double maximumoverlap=qpoverlaps.col(j).maxCoeff(&maxi);
                state_index.push_back(maxi);
                CTP_LOG(ctp::logDEBUG, *_log) << " [" << maxi << " , " << j << "]: " <<maximumoverlap << flush;
              }
            }
          }

          if (_has_osc_filter) {
            // go through list of singlets
            const std::vector<double>oscs = orb_iter_input.Oscillatorstrengths();
            for (unsigned i = 0; i < oscs.size(); i++) {
              if (oscs[i] > _osc_threshold) state_index.push_back(i);
            }
          } else {
            const VectorXfd & energies = (_type == "singlet")
                    ? orb_iter_input.BSESingletEnergies() : orb_iter_input.BSETripletEnergies();
            for (unsigned i = 0; i < energies.size(); i++) {
              state_index.push_back(i);
            }
          }


          // filter according to charge transfer, go through list of excitations in _state_index
          if (_has_dQ_filter) {
            std::vector<int> state_index_copy;
            const std::vector< Eigen::VectorXd >& dQ_frag = (_type == "singlet")
                    ? orb_iter_input.getFragmentChargesSingEXC() : orb_iter_input.getFragmentChargesTripEXC();
            for (unsigned i = 0; i < state_index.size(); i++) {
              if (std::abs(dQ_frag[state_index[i]](0)) > _dQ_threshold) {
                state_index_copy.push_back(state_index[i]);
              }
            }
            state_index = state_index_copy;
          } else if (_has_loc_filter) {
            std::vector<int> state_index_copy;
            const std::vector< Eigen::VectorXd >& popE = (_type == "singlet")
                    ? orb_iter_input.getFragment_E_localisation_singlet() : orb_iter_input.getFragment_E_localisation_triplet();
            const std::vector< Eigen::VectorXd >& popH = (_type == "singlet")
                    ? orb_iter_input.getFragment_H_localisation_singlet() : orb_iter_input.getFragment_H_localisation_triplet();
            if (_localiseonA) {
              for (unsigned i = 0; i < state_index.size(); i++) {
                if (popE[state_index[i]](0) > _loc_threshold && popH[state_index[i]](0) > _loc_threshold) {
                  state_index_copy.push_back(state_index[i]);
                }
              }
            } else {
              for (unsigned i = 0; i < state_index.size(); i++) {
                if (popE[state_index[i]](1) > _loc_threshold && popH[state_index[i]](1) > _loc_threshold) {
                  state_index_copy.push_back(state_index[i]);
                }
              }
            }
            state_index = state_index_copy;
          }


          if (state_index.size() < 1) {
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " WARNING: FILTER yielded no state. Taking lowest excitation" << flush;
            state_index.push_back(0);
          } else {
            if (_type == "quasiparticle") {
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filter yielded QP index: " << state_index[_state - 1 - orb_iter_input.getGWAmin()] << flush;
            } else {
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filter yielded state" << _type << ":" << state_index[_state - 1] + 1 << flush;
            }
          }
          // - output its energy
          if (_type == "singlet") {
            energy_ex = orb_iter_input.BSESingletEnergies()[state_index[_state - 1]] * tools::conv::hrt2ev; // to eV
          } else if (_type == "triplet") {
            energy_ex = orb_iter_input.BSETripletEnergies()[state_index[_state - 1]] * tools::conv::hrt2ev; // to eV
          } else if (_type == "quasiparticle") {
            if (_state > orb_iter_input.getNumberOfElectrons()) {
              // unoccupied QPs: E_a = E_0 + eps_l
              energy_ex = orb_iter_input.QPdiagEnergies()[state_index[_state - 1 - orb_iter_input.getGWAmin()]] * tools::conv::hrt2ev; // to eV
            } else {
              // occupied QPs: E_c = E_0 - eps_h
              energy_ex = -1.0 * orb_iter_input.QPdiagEnergies()[state_index[_state - 1 - orb_iter_input.getGWAmin()]] * tools::conv::hrt2ev; // to eV
            }
          }

        }

        if (!_static_qmmm) {
          Density2Charges(_state);
        }

      } //_do_gwbse


      /* new ESP fit only required for
       * - polarizable QMMM
       * AND
       * - GWBSE or DFT with internal DFTENGINE
       */
      if (!_static_qmmm && _qmpack->getPackageName() == "xtp" && !_do_gwbse) {
        Density2Charges();
      } // for polarized QMMM

      if (tools::globals::verbose) {
        CTP_LOG(ctp::logDEBUG, *_log) << "Calculated partial charges" << flush;
        for (const QMAtom* atom : orb_iter_input.QMAtoms()) {
          CTP_LOG(ctp::logDEBUG, *_log) << atom->getType() << " " << atom->getPartialcharge() << flush;
        }
      }
      
      // EXTRACT & SAVE QM ENERGIES
      double energy_sf = orb_iter_input.getSelfEnergy();
      double energy_qmsf = orb_iter_input.getQMEnergy();
      double energy_qm = energy_qmsf - energy_sf;
      thisIter->setQMSF(energy_qm, energy_sf, energy_ex);
      _job->setEnergy_QMMM(thisIter->getQMEnergy(), thisIter->getGWBSEEnergy(), thisIter->getSFEnergy(),
              thisIter->getQMMMEnergy());

      // EXTRACT & SAVE QMATOM DATA
      std::vector< QMAtom* > &atoms = orb_iter_input.QMAtoms();

      thisIter->UpdatePosChrgFromQMAtoms(atoms, _job->getPolarTop()->QM0());
      
      if (_do_gdma) {
        RunGDMA(thisIter, runFolder);
      } 

      // Update state variable
      if (_type == "quasiparticle" || _has_overlap_filter) {
        _state = state_index[ _state - 1 - orb_iter_input.getGWAmin() ] + 1 + orb_iter_input.getGWAmin();
      }

      // serialize this iteration
      if (_do_archive) {
        // save orbitals
        std::string orb_file = runFolder + "/system.orb";
        CTP_LOG(ctp::logDEBUG, *_log) << "Archiving data to " << orb_file << flush;
        orb_iter_input.WriteToCpt(orb_file);
      }

      CTP_LOG(ctp::logINFO, *_log)
              << format("Summary - iteration %1$d:") % (iterCnt + 1) << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... QM Size  = %1$d atoms") % int(atoms.size()) << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... E(QM)    = %1$+4.9e") % thisIter->getQMEnergy() << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... E(GWBSE) = %1$+4.9e") % thisIter->getGWBSEEnergy() << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... E(SF)    = %1$+4.9e") % thisIter->getSFEnergy() << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... E(FM)    = %1$+4.9e") % thisIter->getFMEnergy() << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... E(MM)    = %1$+4.9e") % thisIter->getMMEnergy() << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... E(QMMM)  = %1$+4.9e") % thisIter->getQMMMEnergy() << flush;
      if (!_static_qmmm) {
        CTP_LOG(ctp::logINFO, *_log)
                << format("... RMS(dR)  = %1$+4.9e") % thisIter->getRMSdR() << flush;
        CTP_LOG(ctp::logINFO, *_log)
                << format("... RMS(dQ)  = %1$+4.9e") % thisIter->getRMSdQ() << flush;
        CTP_LOG(ctp::logINFO, *_log)
                << format("... SUM(dQ)  = %1$+4.9e") % thisIter->getSUMdQ() << flush;
      }
      // CLEAN DIRECTORY
     

      return 0;
    }


    void QMMachine::Density2Charges(int excitedstate_index) {

      Eigen::MatrixXd DMATGS = orb_iter_input.DensityMatrixGroundState();

      Eigen::MatrixXd DMAT_tot = DMATGS; // Ground state + hole_contribution + electron contribution

      if (excitedstate_index > 0) {
        if (_type == "singlet" && _type == "triplet") {
          std::vector<Eigen::MatrixXd > DMAT = orb_iter_input.DensityMatrixExcitedState(_type, excitedstate_index);
          DMAT_tot = DMAT_tot - DMAT[0] + DMAT[1]; // Ground state + hole_contribution + electron contribution
        } else if (_type == "quasiparticle") {

          Eigen::MatrixXd DMATQP = orb_iter_input.DensityMatrixQuasiParticle(excitedstate_index);
          if (excitedstate_index > orb_iter_input.getHomo()) {
            DMAT_tot = DMAT_tot + DMATQP;
          } else {
            DMAT_tot = DMAT_tot - DMATQP;
          }
        }
      }

      
      int iter = _iters.size() - 1;
      Eigen::MatrixXd DMAT_mixed;
      if (iter == 0) {
        DMAT_mixed = DMAT_tot;
      } else {
        DMAT_mixed = _alpha * _DMAT_old + (1 - _alpha) * DMAT_tot;
      }

      _DMAT_old = DMAT_mixed;
      BasisSet dftbs;
      if (orb_iter_input.getDFTbasis() != "") {
        dftbs.LoadBasisSet(orb_iter_input.getDFTbasis());
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " << orb_iter_input.getDFTbasis() << flush;
      }else{
        throw std::runtime_error("Basisset in orb_iter_input not set");
      }
      // fill DFT AO basis by going through all atoms
      AOBasis dftbasis;
      dftbasis.AOBasisFill(dftbs, orb_iter_input.QMAtoms());
      Espfit esp = Espfit(_log);
      esp.Fit2Density(orb_iter_input.QMAtoms(), DMAT_mixed, dftbasis, "medium");

      return;
    }


    QMMIter *QMMachine::CreateNewIter() {

      QMMIter *newIter = new QMMIter(_iters.size());
      this->_iters.push_back(newIter);
      return newIter;
    }


    bool QMMachine::hasConverged() {

    bool  convg_dR = false;
    bool  convg_dQ = false;
    bool  convg_dE_QM = false;
    bool  convg_dE_MM = false;

      if (_iters.size() > 1) {

        QMMIter *iter_0 = _iters[_iters.size() - 2];
        QMMIter *iter_1 = _iters[_iters.size() - 1];

        double dR = iter_1->getRMSdR();
        double dQ = iter_1->getRMSdQ();
        double dE_QM = iter_1->getQMEnergy() - iter_0->getQMEnergy();
        double dE_MM = iter_1->getMMEnergy() - iter_0->getMMEnergy();

        CTP_LOG(ctp::logINFO, *_log)
                << format("... dE_QM  = %1$+4.9e") % dE_QM << flush;
        CTP_LOG(ctp::logINFO, *_log)
                << format("... dE_MM  = %1$+4.9e") % dE_MM << flush;

        if (dR <= _crit_dR) convg_dR = true;
        if (dQ <= _crit_dQ) convg_dQ = true;
        if (dE_QM * dE_QM <= _crit_dE_QM * _crit_dE_QM) convg_dE_QM = true;
        if (dE_MM * dE_MM <= _crit_dE_MM * _crit_dE_MM) convg_dE_MM = true;
      }

      _isConverged = ((convg_dR && convg_dQ) && (convg_dE_QM && convg_dE_MM));

      CTP_LOG(ctp::logINFO, *_log)
              << format("... Convg dR = %s") % (convg_dR ? "true" : "false") << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... Convg dQ = %s") % (convg_dQ ? "true" : "false") << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... Convg QM = %s") % (convg_dE_QM ? "true" : "false") << flush;
      CTP_LOG(ctp::logINFO, *_log)
              << format("... Convg MM = %s") % (convg_dE_MM ? "true" : "false") << flush;

      return _isConverged;
    }




  }
}
