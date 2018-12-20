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
            Property *opt, string sfx)
    : _job(job), _xind(xind), _qmpack(qmpack), 
    _isConverged(false) {

      string key = sfx + ".qmmmconvg";
      _crit_dR = opt->ifExistsReturnElseReturnDefault<double>(key + ".dR", 0.01); //nm
      _crit_dQ = opt->ifExistsReturnElseReturnDefault<double>(key + ".dQ", 0.01); //e
      _crit_dE_QM = opt->ifExistsReturnElseReturnDefault<double>(key + ".dQdE_QM", 0.001); //eV
      _crit_dE_MM = opt->ifExistsReturnElseReturnDefault<double>(key + ".dE_MM", _crit_dE_QM); //eV
      _maxIter = opt->ifExistsReturnElseReturnDefault<int>(key + ".max_iter", 32);
      _alpha = opt->ifExistsReturnElseReturnDefault<double>(key + ".mixing", 0.0);

      key = sfx;
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
        string _gwbse_xml = opt->get(key + ".gwbse_options").as<string> ();
        load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());
         key = sfx + ".gwbse.filter";
         if(opt->exists(key)){
           tools::Property& filterprop=opt->get(key);
           _filter.Initialize(filterprop);
         }
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
      boost::filesystem::path arg_path;
      string frame_dir = "frame_" + boost::lexical_cast<string>(_job->getTop()->getDatabaseId());
      string jobFolder = "job_" + boost::lexical_cast<string>(_job->getId())
              + "_" + _job->getTag();
      string work_dir = (arg_path / "QMMM" / frame_dir / jobFolder).c_str();
      bool created = boost::filesystem::create_directories(work_dir);
      if (created) {
        CTP_LOG(ctp::logINFO, *_log) << "Created directory " << work_dir << flush;
      }
      
      vector<string> split;
      Tokenizer toker(_job->getTag(), ":");
      toker.ToVector(split);
      _initialstate=QMState(split.back());
      if(_initialstate.Type().isExciton() || _initialstate.Type().isGWState()){
       _filter.setInitialState(_initialstate);
       _filter.setLogger(_log);
       _do_gwbse=true;
      }
      
      _qmpack->setCharge(chrg);
      _qmpack->setSpin(spin);

      int iterCnt = 0;
      int iterMax = _maxIter;
      for (; iterCnt < iterMax; ++iterCnt) {

        // check for polarized QM/MM convergence
        int info = Iterate(work_dir, iterCnt);
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
        throw std::runtime_error(" Invalid QMPackage: " + _qmpack->getPackageName()+ " Gaussian 03 only!");
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
      CTP_LOG(ctp::logDEBUG, *_log) << "Running GWBSE "<< flush;
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
      ctp::Logger dft_logger(ctp::logDEBUG);
      dft_logger.setMultithreading(false);
      dft_logger.setPreface(ctp::logINFO, (format("\nGWBSE INF ...")).str());
      dft_logger.setPreface(ctp::logERROR, (format("\nGWBSE ERR ...")).str());
      dft_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...")).str());
      dft_logger.setPreface(ctp::logDEBUG, (format("\nGWBSE DBG ...")).str());
      _qmpack->setLog(&dft_logger);
      
      _qmpack->setMultipoleBackground(MultipolesBackground);
      
      // setting RUNDIR for the external QMPackages, dummy for internal
      _qmpack->setRunDir(runFolder);
      
      CTP_LOG(ctp::logDEBUG, *_log) << "Writing input file " << runFolder << flush;
      CTP_LOG(ctp::logDEBUG, *_log) << "Running DFT " << flush;
      _qmpack->WriteInputFile(orb_iter_input);
      
      orb_iter_input.WriteXYZ(runFolder + "/system_qm.xyz");
      
      _qmpack->Run();
      
      /* Parse information from the LOGFILE into orbitals, if
       * external QMPackage is run. Dummy for internal DFTENGINE.
       * The QMPackage's MOcoefficients are not automatically
       * parsed in DFT-only calculations and ESP fits for
       * polarized QMMM are expected in the LOGFILE.
       * IF the internal ESPFITs should be used, the MOcoefficients
       * need to be parsed too.
       */
      bool success1= _qmpack->ParseLogFile(orb_iter_input);
   
      bool success2=_qmpack->ParseOrbitalsFile(orb_iter_input);
      _qmpack->CleanUp();
      
      ofstream ofs;
      string dft_logfile = runFolder + "/dft.log";
      ofs.open(dft_logfile.c_str(), ofstream::out);
      if (!ofs.is_open()) {
        throw runtime_error("Bad file handle: " + dft_logfile);
      }
      ofs << dft_logger << endl;
      ofs.close();
      _qmpack->setLog(_log);
        return success1&&success2;
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
       

double energy_ex=0.0;
      if (_do_gwbse) {
        if (_initialstate.Type()==QMStateType::Gstate) {
         CTP_LOG(ctp::logDEBUG, *_log) <<"QMMMachine no GWBSE necessary for " + _initialstate.ToLongString();
        }
        RunGWBSE(runFolder);
        CTP_LOG(ctp::logDEBUG, *_log) << "Excited state via GWBSE: " << flush;
        CTP_LOG(ctp::logDEBUG, *_log) << "state:" << _initialstate.ToLongString() << flush;
        _filter.PrintInfo();
       
        QMState newstate= _filter.CalcStateAndUpdate(orb_iter_input);
        energy_ex=orb_iter_input.getExcitedStateEnergy(newstate)* tools::conv::hrt2ev;    

        if (!_static_qmmm) {
          Density2Charges(newstate);
        }

      } //_do_gwbse


      /* new ESP fit only required for
       * - polarizable QMMM
       * AND
       * - GWBSE or DFT with internal DFTENGINE
       */
      if (!_static_qmmm && _qmpack->getPackageName() == "xtp" && !_do_gwbse) {
        QMState gsstate=QMState(QMStateType::Gstate,0,false);
        Density2Charges(gsstate);
      } // for polarized QMMM

      if (tools::globals::verbose) {
        CTP_LOG(ctp::logDEBUG, *_log) << "Calculated partial charges" << flush;
        for (const QMAtom* atom : orb_iter_input.QMAtoms()) {
          CTP_LOG(ctp::logDEBUG, *_log) << atom->getType() << " " << atom->getPartialcharge() << flush;
        }
      }
      
      // EXTRACT & SAVE QM ENERGIES
      double energy_sf = orb_iter_input.getSelfEnergy()*tools::conv::hrt2ev;
      double energy_qmsf = orb_iter_input.getQMEnergy()*tools::conv::hrt2ev;
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


    void QMMachine::Density2Charges(const QMState& state) {

      Eigen::MatrixXd DMAT = orb_iter_input.DensityMatrixFull(state);
      
      int iter = _iters.size() - 1;
      Eigen::MatrixXd DMAT_mixed;
      if (iter == 0) {
        DMAT_mixed = DMAT;
      } else {
        DMAT_mixed = _alpha * _DMAT_old + (1 - _alpha) * DMAT;
      }

      _DMAT_old = DMAT_mixed;
      BasisSet dftbs;
      if (orb_iter_input.getDFTbasisName() != "") {
        dftbs.LoadBasisSet(orb_iter_input.getDFTbasisName());
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " << orb_iter_input.getDFTbasisName() << flush;
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
