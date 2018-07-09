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

#include <string>

#include "idft.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmpackagefactory.h>


using namespace std;
using boost::format;
using namespace boost::filesystem;


namespace votca {
  namespace xtp {

    // +++++++++++++++++++++++++++++ //
    // IDFT MEMBER FUNCTIONS         //
    // +++++++++++++++++++++++++++++ //

    void IDFT::Initialize(votca::tools::Property* options) {

      _energy_difference = 0.0;

      _do_input = false;
      _do_run = false;
      _do_parse = false;
      _do_project = false;
      _do_trim = false;
      _do_extract = false;

      _store_orbitals = false;
      _store_overlap = false;
      _store_integrals = false;

      // update options with the VOTCASHARE defaults   
      UpdateWithDefaults(options, "xtp");
      ParseOptionsXML(options);

      // register all QM packages (Gaussian, turbomole, etc))
      QMPackageFactory::RegisterAll();

    }

    void IDFT::ParseOptionsXML(tools::Property *options) {

      // Orbitals are in fort.7 file; number of electrons in .log file

      std::string key = "options." + Identify();
      _energy_difference = options->get(key + ".degeneracy").as< double > ();


      std::string _tasks_string = options->get(key + ".tasks").as<std::string> ();
      if (_tasks_string.find("input") != std::string::npos) _do_input = true;
      if (_tasks_string.find("run") != std::string::npos) _do_run = true;
      if (_tasks_string.find("parse") != std::string::npos) _do_parse = true;
      if (_tasks_string.find("project") != std::string::npos) _do_project = true;
      if (_tasks_string.find("trim") != std::string::npos) _do_trim = true;
      if (_tasks_string.find("extract") != std::string::npos) _do_extract = true;

      std::string _store_string = options->get(key + ".store").as<std::string> ();
      if (_store_string.find("orbitals") != std::string::npos) _store_orbitals = true;
      if (_store_string.find("overlap") != std::string::npos) _store_overlap = true;
      if (_store_string.find("integrals") != std::string::npos) _store_integrals = true;

    // read linker groups
    std::string linker = options->ifExistsReturnElseReturnDefault<std::string>(key + ".linker_names", "");
    Tokenizer toker(linker, ",");
    toker.ToVector(_linker_names);
 
    
      _max_occupied_levels = options->get(key + ".levels").as<int> ();
      _max_unoccupied_levels = _max_occupied_levels;

      if (options->exists(key + ".trim")) {
        _trim_factor = options->get(key + ".trim").as< int > ();
      } else {
        _trim_factor = -1;
        if (!_do_trim) {
          cout << "WARNING: Are you sure you do not want to trim your orbitals, enable trimming by adding ""trim"" to tasks, default for trimming is -1 e.g. trimming to HOMO/LUMO respectively " << endl;
          ;
        }
      }


      std::string _package_xml = options->get(key + ".dftpackage").as<std::string> ();
      //cout << endl << "... ... Parsing " << _package_xml << endl ;

      load_property_from_xml(_package_options, _package_xml.c_str());

      key = "package";
      _package = _package_options.get(key + ".name").as<std::string> ();

      key = "options." + Identify();

      if (options->exists(key + ".job_file")) {
        _jobfile = options->get(key + ".job_file").as<std::string>();
      } else {
        throw std::runtime_error("Job-file not set. Abort.");
      }


    }

    void IDFT::LoadOrbitals(std::string file_name, Orbitals* orbitals, xtp::Logger *log) {

      XTP_LOG(xtp::logDEBUG, *log) << "Loading " << file_name << flush;
      try {
        orbitals->ReadFromCpt(file_name);
      } catch(std::runtime_error& error){
        XTP_LOG(xtp::logERROR, *log) << "Failed loading orbitals from " << file_name << flush;
      }
    }

    xtp::Job::JobResult IDFT::EvalJob(xtp::Topology *top, xtp::Job *job, xtp::QMThread *opThread) {

      std::string idft_work_dir = "OR_FILES";
      std::string edft_work_dir = "OR_FILES";
      std::string frame_dir = "frame_" + boost::lexical_cast<std::string>(top->getDatabaseId());

      bool _run_status = false;
      bool _parse_log_status = false;
      bool _parse_orbitals_status = false;
      bool _calculate_integrals = false;
      stringstream sout;
      std::string output;

      int HOMO_A;
      int HOMO_B;
      int LUMO_A;
      int LUMO_B;
      Orbitals _orbitalsA, _orbitalsB;
      DFTcoupling dftcoupling;




      // report back to the progress observer
      xtp::Job::JobResult jres = xtp::Job::JobResult();

      // get the logger from the thread
      xtp::Logger* pLog = opThread->getLogger();

      // get the information about the job executed by the thread
      int _job_ID = job->getId();
      Property _job_input = job->getInput();
      list<Property*> segment_list = _job_input.Select("segment");
      int ID_A = segment_list.front()->getAttribute<int>("id");
      std::string type_A = segment_list.front()->getAttribute<std::string>("type");
      int ID_B = segment_list.back()->getAttribute<int>("id");
      std::string type_B = segment_list.back()->getAttribute<std::string>("type");

      // set the folders 
      std::string _pair_dir = (format("%1%%2%%3%%4%%5%") % "pair" % "_" % ID_A % "_" % ID_B).str();

      path arg_path, arg_pathA, arg_pathB, arg_pathAB;

      std::string orbFileA = (arg_pathA / edft_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_A % ".orb").str()).c_str();
      std::string orbFileB = (arg_pathB / edft_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_B % ".orb").str()).c_str();
      std::string orbFileAB = (arg_pathAB / idft_work_dir / "pairs" / frame_dir / (format("%1%%2%%3%%4%%5%") % "pair_" % ID_A % "_" % ID_B % ".orb").str()).c_str();
      std::string _orb_dir = (arg_path / idft_work_dir / "pairs" / frame_dir).c_str();

      xtp::Segment *seg_A = top->getSegment(ID_A);
      assert(seg_A->getName() == type_A);

      xtp::Segment *seg_B = top->getSegment(ID_B);
      assert(seg_B->getName() == type_B);

      XTP_LOG(xtp::logINFO, *pLog) << xtp::TimeStamp() << " Evaluating pair "
              << _job_ID << " [" << ID_A << ":" << ID_B << "] out of " <<
              (top->NBList()).size() << flush;

      std::string _qmpackage_work_dir = (arg_path / idft_work_dir / _package / frame_dir / _pair_dir).c_str();
      // get the corresponding object from the QMPackageFactory
      QMPackage *_qmpackage = QMPackages().Create(_package);
      // set a log file for the package
      _qmpackage->setLog(pLog);
      // set the run dir 
      _qmpackage->setRunDir(_qmpackage_work_dir);
      // get the package options
      _qmpackage->Initialize(&_package_options);


      // if asked, prepare the input files
      if (_do_input) {
        boost::filesystem::create_directories(_qmpackage_work_dir);
        Orbitals *_orbitalsAB = NULL;
        if (_qmpackage->GuessRequested()) { // do not want to do an SCF loop for a dimer

          if (_qmpackage->getPackageName() == "orca") {
            XTP_LOG(xtp::logINFO, *pLog) << "Copying monomer .gbw files to pair folder" << flush;
            std::string gbwFileA = (arg_pathA / edft_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_A % ".gbw").str()).c_str();
            std::string gbwFileB = (arg_pathB / edft_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_B % ".gbw").str()).c_str();
            std::string gbwFileA_workdir = (path(_qmpackage_work_dir) / "molA.gbw").c_str();
            std::string gbwFileB_workdir = (path(_qmpackage_work_dir) / "molB.gbw").c_str();
            boost::filesystem::copy_file(gbwFileA, gbwFileA_workdir, boost::filesystem::copy_option::overwrite_if_exists);
            boost::filesystem::copy_file(gbwFileB, gbwFileB_workdir, boost::filesystem::copy_option::overwrite_if_exists);

          } else {
            XTP_LOG(xtp::logINFO, *pLog) << "Guess requested, reading molecular orbitals" << flush;
            Orbitals _orbitalsA, _orbitalsB;
            _orbitalsAB = new Orbitals();
            // load the corresponding monomer orbitals and prepare the dimer guess 

            // failed to load; wrap-up and finish current job
           
            try {
              _orbitalsA.ReadFromCpt(orbFileA);
            } catch (std::runtime_error& error) {
              XTP_LOG(xtp::logERROR, *pLog) << "Do input: failed loading orbitals from " << orbFileA << flush;
              cout << *pLog;
              output += "failed on " + orbFileA;
              jres.setOutput(output);
              jres.setStatus(xtp::Job::FAILED);
              delete _qmpackage;
              return jres;
            }

            
            try {
              _orbitalsB.ReadFromCpt(orbFileB);
            } catch (std::runtime_error& error) {
              XTP_LOG(xtp::logERROR, *pLog) << "Do input: failed loading orbitals from " << orbFileB << flush;
              cout << *pLog;
              output += "failed on " + orbFileB;
              jres.setOutput(output);
              jres.setStatus(xtp::Job::FAILED);
              delete _qmpackage;
              return jres;
            }
            XTP_LOG(xtp::logERROR, *pLog) << "Writing guess from monomer orbitals" << flush;
            Orbitals::PrepareGuess(&_orbitalsA, &_orbitalsB, _orbitalsAB);
          }
        }
        // if a pair object is available, take into account PBC, otherwise write as is
        xtp::QMNBList* nblist = &top->NBList();
        xtp::QMPair* pair = nblist->FindPair(seg_A, seg_B);

        if (pair == NULL) {
          vector < xtp::Segment* > segments;
          segments.push_back(seg_A);
          segments.push_back(seg_B);
          XTP_LOG(xtp::logWARNING, *pLog) << "PBCs are not taken into account when writing the coordinate file!" << flush;
          _qmpackage->WriteInputFile(segments, _orbitalsAB);
        } else {
            _qmpackage->WriteInputFilePBC(pair, _orbitalsAB, _linker_names);
        }

        delete _orbitalsAB;
      } // end of the input


      // run the executable
      if (_do_run) {
        _run_status = _qmpackage->Run();
        if (!_run_status) {
          output += "run failed; ";
          XTP_LOG(xtp::logERROR, *pLog) << _qmpackage->getPackageName() << " run failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(xtp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }
      } // end of the run


      // This will be later used to write orbitals of the dimer to a file 
      // SOMETHING TO CLEANUP
      Orbitals _orbitalsAB;
      // parse the log/orbitals files
      if (_do_parse) {
        _parse_log_status = _qmpackage->ParseLogFile(&_orbitalsAB);

        if (!_parse_log_status) {
          output += "log incomplete; ";
          XTP_LOG(xtp::logERROR, *pLog) << "LOG parsing failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(xtp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }



        _parse_orbitals_status = _qmpackage->ParseOrbitalsFile(&_orbitalsAB);

        if (!_parse_orbitals_status) {
          output += "Orbitals parsing failed; ";
          XTP_LOG(xtp::logERROR, *pLog) << "Orbitals parsing failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(xtp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }
      } // end of the parse orbitals/log

      // orbital file used to archive parsed data
      std::string _pair_file = (format("%1%%2%%3%%4%%5%") % "pair_" % ID_A % "_" % ID_B % ".orb").str();
      Eigen::MatrixXd _JAB;


      Property _job_summary;

      // Orbitals _orbitalsA, _orbitalsB;

      int _degAH = 1;
      int _degAL = 1;
      int _degBH = 1;
      int _degBL = 1;
      if (_do_project) {

        // orbitals must be loaded from a file
        if (!_do_parse) LoadOrbitals(orbFileAB, &_orbitalsAB, pLog);


        // failed to load; wrap-up and finish current job
        try {
          _orbitalsA.ReadFromCpt(orbFileA);
        } catch (std::runtime_error& error) {
          XTP_LOG(xtp::logERROR, *pLog) << "Failed loading orbitals from " << orbFileA << flush;
          cout << *pLog;
          output += "failed on " + orbFileA;
          jres.setOutput(output);
          jres.setStatus(xtp::Job::FAILED);
          delete _qmpackage;
          return jres;
        } 
        try {
          _orbitalsB.ReadFromCpt(orbFileB);
        } catch (std::runtime_error& error) {
          XTP_LOG(xtp::logERROR, *pLog) << "Failed loading orbitals from " << orbFileB << flush;
          cout << *pLog;
          output += "failed on " + orbFileB;
          jres.setOutput(output);
          jres.setStatus(xtp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }



        if (_trim_factor == -1) {

          // find degeneracy of HOMOs and LUMOs
          std::vector<int> list_levelsAH = (*_orbitalsA.getDegeneracy(_orbitalsA.getNumberOfElectrons() - 1, _energy_difference));
          _degAH = list_levelsAH.size();
          std::vector<int> list_levelsAL = (*_orbitalsA.getDegeneracy(_orbitalsA.getNumberOfElectrons(), _energy_difference));
          _degAL = list_levelsAL.size();

          std::vector<int> list_levelsBH = (*_orbitalsB.getDegeneracy(_orbitalsB.getNumberOfElectrons() - 1, _energy_difference));
          _degBH = list_levelsBH.size();
          std::vector<int> list_levelsBL = (*_orbitalsB.getDegeneracy(_orbitalsB.getNumberOfElectrons(), _energy_difference));
          _degBL = list_levelsBL.size();
        }

        if (_do_trim) {


          if (_trim_factor == -1) {
            XTP_LOG(xtp::logDEBUG, *pLog) << "Trimming orbitals to minimal set " << flush;
            XTP_LOG(xtp::logDEBUG, *pLog) << "HOMO(A)-" << _degAH << " to " << "LUMO(A)+" << _degAL << flush;
            _orbitalsA.Trim(_degAH, _degAL);
            XTP_LOG(xtp::logDEBUG, *pLog) << "HOMO(B)-" << _degBH << " to " << "LUMO(B)+" << _degBL << flush;
            _orbitalsB.Trim(_degBH, _degBL);

          } else {


            XTP_LOG(xtp::logDEBUG, *pLog) << "Trimming virtual orbitals A:"
                    << _orbitalsA.getNumberOfLevels() - _orbitalsA.getNumberOfElectrons() << "->"
                    << _orbitalsA.getNumberOfElectrons()*(_trim_factor - 1);
            _orbitalsA.Trim(_trim_factor);

            XTP_LOG(xtp::logDEBUG, *pLog) << " B:"
                    << _orbitalsB.getNumberOfLevels() - _orbitalsB.getNumberOfElectrons() << "->"
                    << _orbitalsB.getNumberOfElectrons()*(_trim_factor - 1) << flush;
            _orbitalsB.Trim(_trim_factor);

          }
        } // _do_trim

        dftcoupling.setLogger(pLog);

        
        _calculate_integrals = dftcoupling.CalculateIntegrals(&_orbitalsA, &_orbitalsB, &_orbitalsAB, &_JAB);

        if (!_calculate_integrals) {
          output += "integrals failed; ";
          XTP_LOG(xtp::logERROR, *pLog) << "Calculating integrals failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(xtp::Job::FAILED);
          return jres;
        }

        HOMO_A = _orbitalsA.getNumberOfElectrons();
        HOMO_B = _orbitalsB.getNumberOfElectrons();
        LUMO_A = HOMO_A + 1;
        LUMO_B = HOMO_B + 1;

        double J_h;
        double J_e;

        if (_trim_factor == -1) {
          J_h = dftcoupling.getCouplingElement(_degAH, _degBH, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference);
          J_e = dftcoupling.getCouplingElement(_degAH + 1, _degBH + 1, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference);
        } else {
          J_h = dftcoupling.getCouplingElement(HOMO_A, HOMO_B, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference);
          J_e = dftcoupling.getCouplingElement(LUMO_A, LUMO_B, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference);
        }
        XTP_LOG(xtp::logINFO, *pLog) << "Couplings h/e " << ID_A << ":" << ID_B << " " << J_h << ":" << J_e << flush;

        // Output the thread run summary and clean the Logger
        XTP_LOG(xtp::logINFO, *pLog) << xtp::TimeStamp() << " Finished evaluating pair " << ID_A << ":" << ID_B << flush;
        //cout << *pLog;


        if (!(_store_orbitals && _do_parse && _parse_orbitals_status)) {
          _store_orbitals = false;
          XTP_LOG(xtp::logINFO, *pLog) << "Not storing orbitals" << flush;
        }
        if (!(_store_overlap && _do_parse && _parse_log_status)) {
          _store_overlap = false;
          XTP_LOG(xtp::logINFO, *pLog) << "Not storing overlap" << flush;
        }
        if (!(_store_integrals && _do_project && _calculate_integrals)) {
          _store_integrals = false;
          XTP_LOG(xtp::logINFO, *pLog) << "Not storing integrals" << flush;
        } else {
          // _orbitalsAB.setIntegrals( &_JAB );
          Eigen::MatrixXd& _JAB_store = _orbitalsAB.MOCouplings();
          _JAB_store = _JAB;
        }
        // save orbitals 
        boost::filesystem::create_directories(_orb_dir);

        XTP_LOG(xtp::logDEBUG, *pLog) << "Saving orbitals to " << _pair_file << flush;

        _orbitalsAB.WriteToCpt(orbFileAB);

      } // end of the projection loop


      if (_do_extract) {
        LoadOrbitals(orbFileAB, &_orbitalsAB, pLog);
        LoadOrbitals(orbFileA, &_orbitalsA, pLog);
        LoadOrbitals(orbFileB, &_orbitalsB, pLog);
        if (_trim_factor == -1) {

          // find degeneracy of HOMOs and LUMOs
          std::vector<int> list_levelsAH = (*_orbitalsA.getDegeneracy(_orbitalsA.getNumberOfElectrons() - 1, _energy_difference));
          _degAH = list_levelsAH.size();
          std::vector<int> list_levelsAL = (*_orbitalsA.getDegeneracy(_orbitalsA.getNumberOfElectrons(), _energy_difference));
          _degAL = list_levelsAL.size();

          std::vector<int> list_levelsBH = (*_orbitalsB.getDegeneracy(_orbitalsB.getNumberOfElectrons() - 1, _energy_difference));
          _degBH = list_levelsBH.size();
          std::vector<int> list_levelsBL = (*_orbitalsB.getDegeneracy(_orbitalsB.getNumberOfElectrons(), _energy_difference));
          _degBL = list_levelsBL.size();

          _orbitalsA.Trim(_degAH, _degAL);
          _orbitalsB.Trim(_degBH, _degBL);

        } else {
          _orbitalsA.Trim(_trim_factor);
          _orbitalsB.Trim(_trim_factor);
        }
        _JAB = _orbitalsAB.MOCouplings();
        HOMO_A = _orbitalsA.getNumberOfElectrons();
        HOMO_B = _orbitalsB.getNumberOfElectrons();
        LUMO_A = HOMO_A + 1;
        LUMO_B = HOMO_B + 1;


      }

      Property *_job_output = &_job_summary.add("output", "");
      Property *_pair_summary = &_job_output->add("pair", "");
      std::string nameA = seg_A->getName();
      std::string nameB = seg_B->getName();
      _pair_summary->setAttribute("idA", ID_A);
      _pair_summary->setAttribute("idB", ID_B);
      _pair_summary->setAttribute("typeA", nameA);
      _pair_summary->setAttribute("typeB", nameB);

      if (_do_project || _do_extract) {
        _pair_summary->setAttribute("homoA", HOMO_A);
        _pair_summary->setAttribute("homoB", HOMO_B);

        if (_trim_factor == -1) {

          // HOMO-HOMO coupling
          double JAB = dftcoupling.getCouplingElement(_degAH, _degBH, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference);
          Property *_overlap_summary = &_pair_summary->add("overlap", boost::lexical_cast<std::string>(JAB));
          double energyA = _orbitalsA.getEnergy(_degAH);
          double energyB = _orbitalsB.getEnergy(_degBH);
          _overlap_summary->setAttribute("orbA", HOMO_A);
          _overlap_summary->setAttribute("orbB", HOMO_B);
          _overlap_summary->setAttribute("jAB", JAB);
          _overlap_summary->setAttribute("eA", energyA);
          _overlap_summary->setAttribute("eB", energyB);

          // LUMO-LUMO coupling
          JAB = dftcoupling.getCouplingElement(_degAH + 1, _degBH + 1, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference);
          _overlap_summary = &_pair_summary->add("overlap", boost::lexical_cast<std::string>(JAB));
          energyA = _orbitalsA.getEnergy(_degAH + 1);
          energyB = _orbitalsB.getEnergy(_degBH + 1);
          _overlap_summary->setAttribute("orbA", LUMO_A);
          _overlap_summary->setAttribute("orbB", LUMO_B);
          _overlap_summary->setAttribute("jAB", JAB);
          _overlap_summary->setAttribute("eA", energyA);
          _overlap_summary->setAttribute("eB", energyB);

        } else {


          for (int levelA = HOMO_A - _max_occupied_levels + 1; levelA <= LUMO_A + _max_unoccupied_levels - 1; ++levelA) {
            for (int levelB = HOMO_B - _max_occupied_levels + 1; levelB <= LUMO_B + _max_unoccupied_levels - 1; ++levelB) {
              Property *_overlap_summary = &_pair_summary->add("overlap", "");
              double JAB = dftcoupling.getCouplingElement(levelA, levelB, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference);
              double energyA = _orbitalsA.getEnergy(levelA);
              double energyB = _orbitalsB.getEnergy(levelB);
              _overlap_summary->setAttribute("orbA", levelA);
              _overlap_summary->setAttribute("orbB", levelB);
              _overlap_summary->setAttribute("jAB", JAB);
              _overlap_summary->setAttribute("eA", energyA);
              _overlap_summary->setAttribute("eB", energyB);
            }
          }
        }

        votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
        sout << iomXML << _job_summary;
        // Flo: moved next line from below
        jres.setOutput(_job_summary);
      }


      // cleanup whatever is not needed
      _qmpackage->CleanUp();
      delete _qmpackage;
      jres.setOutput(_job_summary);
      jres.setStatus(xtp::Job::COMPLETE);

      return jres;
    }

    void IDFT::WriteJobFile(xtp::Topology *top) {

      cout << endl << "... ... Writing job file " << flush;
      std::ofstream ofs;
      ofs.open(_jobfile.c_str(), std::ofstream::out);
      if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);


      xtp::QMNBList::iterator pit;
      xtp::QMNBList &nblist = top->NBList();

      int jobCount = 0;
      if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
      }

      ofs << "<jobs>" << endl;
      std::string tag = "";

      for (pit = nblist.begin(); pit != nblist.end(); ++pit) {

        int id1 = (*pit)->Seg1()->getId();
        std::string name1 = (*pit)->Seg1()->getName();
        int id2 = (*pit)->Seg2()->getId();
        std::string name2 = (*pit)->Seg2()->getName();

        int id = ++jobCount;

        Property Input;
        Property *pInput = &Input.add("input", "");
        Property *pSegment = &pInput->add("segment", boost::lexical_cast<std::string>(id1));
        pSegment->setAttribute<string>("type", name1);
        pSegment->setAttribute<int>("id", id1);

        pSegment = &pInput->add("segment", boost::lexical_cast<string>(id2));
        pSegment->setAttribute<string>("type", name2);
        pSegment->setAttribute<int>("id", id2);

        xtp::Job job(id, tag, Input, xtp::Job::AVAILABLE);
        job.ToStream(ofs, "xml");

      }

      // CLOSE STREAM
      ofs << "</jobs>" << endl;
      ofs.close();

      cout << endl << "... ... In total " << jobCount << " jobs" << flush;

    }

    /** 
     * Imports electronic couplings with superexchange
     */

    void IDFT::ReadJobFile(xtp::Topology *top) {

      Property xml;

      vector<Property*> records;

      // gets the neighborlist from the topology
      xtp::QMNBList &nblist = top->NBList();
      int _number_of_pairs = nblist.size();
      int _current_pairs = 0;
      int _incomplete_jobs = 0;

      // output using logger
      xtp::Logger _log;
      _log.setReportLevel(xtp::logINFO);

      // generate lists of bridges for superexchange pairs
      nblist.GenerateSuperExchange();

      // load the QC results in a vector indexed by the pair ID
      load_property_from_xml(xml, _jobfile);
      list<Property*> jobProps = xml.Select("jobs.job");

      records.resize(jobProps.size() + 1);

      // loop over all jobs = pair records in the job file
      for (list<Property*> ::iterator it = jobProps.begin(); it != jobProps.end(); ++it) {

        // if job produced an output, then continue with analysis
        if ((*it)->exists("output") && (*it)->exists("output.pair")) {

          // get the output records
          Property poutput = (*it)->get("output.pair");
          // id's of two segments of a pair
          int idA = poutput.getAttribute<int>("idA");
          int idB = poutput.getAttribute<int>("idB");
          // segments which correspond to these ids           
          xtp::Segment *segA = top->getSegment(idA);
          xtp::Segment *segB = top->getSegment(idB);
          // pair that corresponds to the two segments
          xtp::QMPair *qmp = nblist.FindPair(segA, segB);

          if (qmp == NULL) { // there is no pair in the neighbor list with this name
            XTP_LOG_SAVE(xtp::logINFO, _log) << "No pair " << idA << ":" << idB << " found in the neighbor list. Ignoring" << flush;
          } else {
            //XTP_LOG(logINFO, _log) << "Store in record: " <<  idA << ":" << idB << flush; 
            records[qmp->getId()] = &((*it)->get("output.pair"));
          }
        }
      } // finished loading from the file


      // loop over all pairs in the neighbor list
      std::cout << "Neighborlist size " << top->NBList().size() << std::endl;
      for (xtp::QMNBList::iterator ipair = top->NBList().begin(); ipair != top->NBList().end(); ++ipair) {

        xtp::QMPair *pair = *ipair;
        xtp::Segment* segmentA = pair->Seg1();
        xtp::Segment* segmentB = pair->Seg2();

        double Jhop2_homo = 0;
        double Jeff_homo = 0;
        double Jeff_lumo = 0;

        cout << "\nProcessing pair " << segmentA->getId() << ":" << segmentB->getId() << flush;

        xtp::QMPair::PairType _ptype = pair->getType();
        Property* pair_property = records[ pair->getId() ];

        if (pair_property) {
          int homoA = pair_property->getAttribute<int>("homoA");
          int homoB = pair_property->getAttribute<int>("homoB");

          // If a pair is of a direct type 
          if (_ptype == xtp::QMPair::Hopping || _ptype == xtp::QMPair::SuperExchangeAndHopping) {
            cout << ":hopping";
            list<Property*> pOverlap = pair_property->Select("overlap");

            for (list<Property*> ::iterator itOverlap = pOverlap.begin(); itOverlap != pOverlap.end(); ++itOverlap) {

              double overlapAB = (*itOverlap)->getAttribute<double>("jAB");
              int orbA = (*itOverlap)->getAttribute<double>("orbA");
              int orbB = (*itOverlap)->getAttribute<double>("orbB");
              // cout << " orbA:orbB " << orbA << ":" << orbB << flush;
              if (orbA == homoA && orbB == homoB) {
                Jeff_homo += overlapAB;
                Jhop2_homo += overlapAB*overlapAB;
              }

              if (orbA == homoA + 1 && orbB == homoB + 1) {
                Jeff_lumo += overlapAB;
              }

            }

          }


          double Jeff2_homo = Jeff_homo*Jeff_homo;
          double Jeff2_lumo = Jeff_lumo*Jeff_lumo;


          cout << " Jhop2_HOMO: " << Jhop2_homo << endl;
          cout << " Jeff2_HOMO: " << Jeff2_homo << " (+" << (Jeff2_homo - Jhop2_homo) / Jhop2_homo * 100 << " %)" << endl;
          // cout << " Jeff2_LUMO: " << Jeff2_lumo << endl;

          cout << endl;

          pair->setJeff2(Jeff2_homo, 1);
          pair->setIsPathCarrier(true, 1);

          pair->setJeff2(Jeff2_lumo, -1);
          pair->setIsPathCarrier(true, -1);

        } else {
          cout << " Job incomplete. Setting couplings to 0.0." << endl;
          _incomplete_jobs += 1;

          pair->setJeff2(0.0, 1);
          //pair->setIsPathCarrier(true, 1);

          pair->setJeff2(0.0, -1);
          //pair->setIsPathCarrier(true, -1);
        }

      }

      XTP_LOG_SAVE(xtp::logINFO, _log) << "Pairs [total:updated] " << _number_of_pairs << ":" << _current_pairs << " Incomplete jobs: " << _incomplete_jobs << flush;
      cout << _log;
      return;
    }



  }
};
