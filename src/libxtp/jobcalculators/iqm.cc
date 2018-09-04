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


#include "iqm.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/split.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/constants.h>
#include <votca/xtp/qminterface.h>
#include <votca/xtp/qmpackagefactory.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
  namespace xtp {

    void IQM::Initialize(votca::tools::Property* options) {


      _do_dft_input = false;
      _do_dft_run = false;
      _do_dft_parse = false;
      _do_dftcoupling = false;
      _do_gwbse = false;
      _do_bsecoupling = false;

      _store_dft = false;
      _store_singlets = false;
      _store_triplets = false;
      _store_ehint = false;
      ParseOptionsXML(*options);

      // register all QM packages (Gaussian, turbomole, etc))
      QMPackageFactory::RegisterAll();
      return;
    }

    void IQM::ParseOptionsXML(votca::tools::Property &opt) {


      // parsing general ibse options
      string key = "options." + Identify();
      // _energy_difference = opt.get( key + ".degeneracy" ).as< double > ();


      // job tasks
      string tasks_string = opt.get(key + ".tasks").as<string> ();
      if (tasks_string.find("input") != std::string::npos) _do_dft_input = true;
      if (tasks_string.find("dft") != std::string::npos) _do_dft_run = true;
      if (tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
      if (tasks_string.find("dftcoupling") != std::string::npos) _do_dftcoupling = true;
      if (tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;
      if (tasks_string.find("bsecoupling") != std::string::npos) _do_bsecoupling = true;

      // storage options
      string store_string = opt.get(key + ".store").as<string> ();
      if (store_string.find("dft") != std::string::npos) _store_dft = true;
      if (store_string.find("singlets") != std::string::npos) _store_singlets = true;
      if (store_string.find("triplets") != std::string::npos) _store_triplets = true;
      if (store_string.find("ehint") != std::string::npos) _store_ehint = true;


      if (_do_dft_input || _do_dft_run || _do_dft_parse) {
        string _package_xml = opt.get(key + ".dftpackage").as<string> ();
        load_property_from_xml(_dftpackage_options, _package_xml.c_str());
        string dftname = "package.name";
        _package = _dftpackage_options.get(dftname).as<string> ();
      }

      // read linker groups
      string linker = opt.ifExistsReturnElseReturnDefault<string>(key + ".linker_names", "");
      Tokenizer toker(linker, ",");
      toker.ToVector(_linker_names);

      if (_do_dftcoupling) {
        _dftcoupling_options = opt.get(key + ".dftcoupling_options");
      }

      if (_do_gwbse) {
        string _gwbse_xml = opt.get(key + ".gwbse_options").as<string> ();
        load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());
      }
      if (_do_bsecoupling) {
        string _coupling_xml = opt.get(key + ".bsecoupling_options").as<string>();
        load_property_from_xml(_bsecoupling_options, _coupling_xml.c_str());
      }


      //options for parsing data into sql file   
      key = "options." + Identify() + ".readjobfile";
      if (opt.exists(key + ".singlet")) {
        string parse_string_s = opt.get(key + ".singlet").as<string> ();
        _singlet_levels = FillParseMaps(parse_string_s);
      }
      if (opt.exists(key + ".triplet")) {
        string parse_string_t = opt.get(key + ".triplet").as<string> ();
        _triplet_levels = FillParseMaps(parse_string_t);
      }

      if (opt.exists(key + ".hole")) {
        string parse_string_h = opt.get(key + ".hole").as<string> ();
        _hole_levels = FillParseMaps(parse_string_h);
      }
      if (opt.exists(key + ".electron")) {
        string parse_string_e = opt.get(key + ".electron").as<string> ();
        _electron_levels = FillParseMaps(parse_string_e);
      }

      // job file specification
      key = "options." + Identify();

      if (opt.exists(key + ".job_file")) {
        _jobfile = opt.get(key + ".job_file").as<string>();
      } else {
        throw std::runtime_error("Job-file not set. Abort.");
      }

      return;
    }

    std::map<std::string, QMState> IQM::FillParseMaps(const string& Mapstring) {
      Tokenizer split_options(Mapstring, ", \t\n");
      std::map<std::string, QMState> type2level;
      for (const string& substring : split_options) {
        std::vector<string> segmentpnumber;
        Tokenizer tok(substring, ":");
        tok.ToVector(segmentpnumber);
        if (segmentpnumber.size() != 2) {
          throw runtime_error("Parser iqm: Segment and exciton labels:" + substring + "are not separated properly");
        }
        QMState state=QMState(segmentpnumber[1]);
        string segmentname = segmentpnumber[0];
        type2level[segmentname] = state;
      }
      return type2level;
    }

    void IQM::addLinkers(std::vector< ctp::Segment* > &segments, ctp::Topology *top) {
      ctp::Segment* seg1 = segments[0];
      ctp::Segment* seg2 = segments[1];
      int moleculeIdSeg1 = seg1->getMolecule()-> getId();
      int moleculeIdSeg2 = seg2->getMolecule()-> getId();
      if (moleculeIdSeg1 == moleculeIdSeg2) {// Check that both segments belong to the same molecule
        int idSeg1 = seg1->getId();
        int idSeg2 = seg2->getId();
        std::vector<ctp::Segment*> segmentsInMolecule = top->getMolecule(moleculeIdSeg1)->Segments();
        for (ctp::Segment* segment : segmentsInMolecule) {
          int idIterator = segment->getId();
          if (idIterator != idSeg1 && idIterator != idSeg2 && isLinker(segment->getName())) {
            segments.push_back(segment);
          }
        }
      }
      return;
    }

    bool IQM::isLinker(const std::string& name) {
      return (std::find(_linker_names.begin(), _linker_names.end(), name) != _linker_names.end());
    }

    void IQM::WriteCoordinatesToOrbitalsPBC(ctp::QMPair& pair, Orbitals& orbitals) {
      ctp::Segment* seg1 = pair.Seg1();
      ctp::Segment* seg2 = pair.Seg2();
      ctp::Segment* ghost = NULL;
      ctp::Topology* _top = seg1->getTopology();
      ctp::vec r1 = seg1->getPos();
      ctp::vec r2 = seg2->getPos();
      ctp::vec _R = _top->PbShortestConnect(r1, r2); // => _R points from 1 to 2

      // Check whether pair formed across periodic boundary
      if (abs(r2 - r1 - _R) > 1e-8) {
        ghost = new ctp::Segment(seg2);
        //ghost->TranslateBy(r1 - r2 + _R); // DO NOT USE THIS METHOD !
        for (ctp::Atom* atom : ghost->Atoms()) {
          atom->setQMPos(atom->getQMPos() + r1 - r2 + _R);
        }
      }
      std::vector< ctp::Segment* > segments;
      segments.push_back(seg1);
      if (ghost) {
        segments.push_back(ghost);
      } else {
        segments.push_back(seg2);
      }
      QMInterface interface;
      orbitals.QMAtoms() = interface.Convert(segments);
      delete ghost;
    }

    void IQM::SetJobToFailed(ctp::Job::JobResult& jres, ctp::Logger* pLog, const string& errormessage) {
      CTP_LOG(ctp::logERROR, *pLog) << errormessage << flush;
      cout << *pLog;
      jres.setError(errormessage);
      jres.setStatus(ctp::Job::FAILED);
    }
    
    void IQM::WriteLoggerToFile(const string& logfile, ctp::Logger& logger){
      std::ofstream ofs;
      ofs.open(logfile.c_str(), std::ofstream::out);
      if (!ofs.is_open()) {
        throw runtime_error("Bad file handle: " + logfile);
      }
      ofs << logger << endl;
      ofs.close();
    }

    ctp::Job::JobResult IQM::EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *opThread) {

      // report back to the progress observer
      ctp::Job::JobResult jres = ctp::Job::JobResult();

      string iqm_work_dir = "OR_FILES";
      string eqm_work_dir = "OR_FILES";
      string frame_dir = "frame_" + boost::lexical_cast<string>(top->getDatabaseId());

      ctp::Logger* pLog = opThread->getLogger();

      // get the information about the job executed by the thread
      int job_ID = job->getId();
      Property job_input = job->getInput();
      list<Property*> segment_list = job_input.Select("segment");
      int ID_A = segment_list.front()->getAttribute<int>("id");
      string type_A = segment_list.front()->getAttribute<string>("type");
      int ID_B = segment_list.back()->getAttribute<int>("id");
      string type_B = segment_list.back()->getAttribute<string>("type");

      // set the folders 
      string pair_dir = (format("%1%%2%%3%%4%%5%") % "pair" % "_" % ID_A % "_" % ID_B).str();
       
      boost::filesystem::path arg_path, arg_pathA, arg_pathB, arg_pathAB;
      
      string orbFileA = (arg_pathA / eqm_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_A % ".orb").str()).c_str();
      string orbFileB = (arg_pathB / eqm_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_B % ".orb").str()).c_str();
      string orbFileAB = (arg_pathAB / iqm_work_dir / "pairs_iqm" / frame_dir / (format("%1%%2%%3%%4%%5%") % "pair_" % ID_A % "_" % ID_B % ".orb").str()).c_str();
      string orb_dir = (arg_path / iqm_work_dir / "pairs_iqm" / frame_dir).c_str();

      ctp::Segment *seg_A = top->getSegment(ID_A);
      ctp::Segment *seg_B = top->getSegment(ID_B);
      ctp::QMNBList* nblist = &top->NBList();
      ctp::QMPair* pair = nblist->FindPair(seg_A, seg_B);

      CTP_LOG(ctp::logINFO, *pLog) << ctp::TimeStamp() << " Evaluating pair "
              << job_ID << " [" << ID_A << ":" << ID_B << "] out of " <<
              (top->NBList()).size() << flush;

      string package_append = _package + "_" + Identify();
      std::vector< ctp::Segment* > segments;
      segments.push_back(seg_A);
      segments.push_back(seg_B);
      string work_dir = (arg_path / iqm_work_dir / package_append / frame_dir / pair_dir).c_str();

      if (_linker_names.size() > 0) {
        addLinkers(segments, top);
      }
      Orbitals orbitalsAB;
      // if a pair object is available and is not linked take into account PBC, otherwise write as is
      if (pair == NULL || segments.size() > 2) {
        if (pair == NULL) {
          CTP_LOG(ctp::logWARNING, *pLog) << "PBCs are not taken into account when writing the coordinate file!" << flush;
        }
        QMInterface interface;
        orbitalsAB.QMAtoms() = interface.Convert(segments);
      } else {
        WriteCoordinatesToOrbitalsPBC(*pair, orbitalsAB);
      }

      if (_do_dft_input || _do_dft_run || _do_dft_parse) {
        string qmpackage_work_dir = (arg_path / iqm_work_dir / package_append / frame_dir / pair_dir).c_str();
        
        ctp::Logger dft_logger(ctp::logDEBUG);
        dft_logger.setMultithreading(false);
        dft_logger.setPreface(ctp::logINFO, (format("\nDFT INF ...")).str());
        dft_logger.setPreface(ctp::logERROR, (format("\nDFT ERR ...")).str());
        dft_logger.setPreface(ctp::logWARNING, (format("\nDFT WAR ...")).str());
        dft_logger.setPreface(ctp::logDEBUG, (format("\nDFT DBG ...")).str());

        QMPackage *qmpackage = QMPackages().Create(_package);
        qmpackage->setLog(&dft_logger);
        qmpackage->setRunDir(qmpackage_work_dir);
        qmpackage->Initialize(_dftpackage_options);


        // if asked, prepare the input files
        if (_do_dft_input) {
          boost::filesystem::create_directories(qmpackage_work_dir);
          if (qmpackage->GuessRequested()) {
            if (_linker_names.size() > 0) {
              throw std::runtime_error("Error: You are using a linker and want "
                      "to use a monomer guess for the dimer. These are mutually exclusive.");
            }

            CTP_LOG(ctp::logINFO, *pLog) << "Guess requested, reading molecular orbitals" << flush;

            if (qmpackage->getPackageName() == "orca") {
              CTP_LOG(ctp::logINFO, *pLog) << "Copying monomer .gbw files to pair folder" << flush;
              string gbwFileA = (arg_pathA / eqm_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_A % ".gbw").str()).c_str();
              string gbwFileB = (arg_pathB / eqm_work_dir / "molecules" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_B % ".gbw").str()).c_str();
              string gbwFileA_workdir = (qmpackage_work_dir / "molA.gbw").c_str();
              string gbwFileB_workdir = (qmpackage_work_dir / "molB.gbw").c_str();
              boost::filesystem::copy_file(gbwFileA, gbwFileA_workdir, boost::filesystem::copy_option::overwrite_if_exists);
              boost::filesystem::copy_file(gbwFileB, gbwFileB_workdir, boost::filesystem::copy_option::overwrite_if_exists);
            } else {
              Orbitals orbitalsB;
              Orbitals orbitalsA;

              try {
                orbitalsA.ReadFromCpt(orbFileA);
              } catch (std::runtime_error& error) {
                SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileA);
                delete qmpackage;
                return jres;
              }

              try {
                orbitalsB.ReadFromCpt(orbFileB);
              } catch (std::runtime_error& error) {
                SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileB);
                delete qmpackage;
                return jres;
              }
              CTP_LOG(ctp::logDEBUG, *pLog) << "Constructing the guess for dimer orbitals" << flush;
              orbitalsAB.PrepareDimerGuess(orbitalsA, orbitalsB);
            }
          } else {
            CTP_LOG(ctp::logINFO, *pLog) << "No Guess requested, starting from DFT starting Guess" << flush;
          }
          qmpackage->WriteInputFile(orbitalsAB);
        }

        if (_do_dft_run) {
          CTP_LOG(ctp::logDEBUG, *pLog) << "Running DFT" << flush;
          bool _run_dft_status = qmpackage->Run(orbitalsAB);
          if (!_run_dft_status) {
            SetJobToFailed(jres, pLog, qmpackage->getPackageName() + " run failed");
            delete qmpackage;
            return jres;
          }
        }

        if (_do_dft_parse) {
          bool parse_log_status = qmpackage->ParseLogFile(orbitalsAB);
          if (!parse_log_status) {
            SetJobToFailed(jres, pLog, "LOG parsing failed");
            delete qmpackage;
            return jres;
          }

          bool parse_orbitals_status = qmpackage->ParseOrbitalsFile(orbitalsAB);

          if (!parse_orbitals_status) {
            SetJobToFailed(jres, pLog, "Orbitals parsing failed");
            delete qmpackage;
            return jres;
          }

        }// end of the parse orbitals/log
        qmpackage->CleanUp();
        delete qmpackage;
        WriteLoggerToFile(work_dir + "/dft.log",dft_logger);
      } else {
        try {
          orbitalsAB.ReadFromCpt(orbFileAB);
        } catch (std::runtime_error& error) {
          SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileAB);
          return jres;
        }
      }
      Property _job_summary;
      Property &job_output = _job_summary.add("output", "");
      if (_do_dftcoupling) {
        DFTcoupling dftcoupling;
        dftcoupling.setLogger(pLog);
        dftcoupling.Initialize(_dftcoupling_options);
        Orbitals orbitalsB;
        Orbitals orbitalsA;

        try {
          orbitalsA.ReadFromCpt(orbFileA);
        } catch (std::runtime_error& error) {
          SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileA);
          return jres;
        }

        try {
          orbitalsB.ReadFromCpt(orbFileB);
        } catch (std::runtime_error& error) {
          SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileB);
          return jres;
        }

        try {
          dftcoupling.CalculateCouplings(orbitalsA, orbitalsB, orbitalsAB);
          dftcoupling.Addoutput(job_output, orbitalsA, orbitalsB);
        } catch (std::runtime_error& error) {
          std::string errormessage(error.what());
          SetJobToFailed(jres, pLog, errormessage);
          return jres;
        }
      }


      // do excited states calculation
      if (_do_gwbse) {
        try {
          CTP_LOG(ctp::logDEBUG, *pLog) << "Running GWBSE" << flush;
          ctp::Logger gwbse_logger(ctp::logDEBUG);
          gwbse_logger.setMultithreading(false);
          gwbse_logger.setPreface(ctp::logINFO, (format("\nGWBSE INF ...")).str());
          gwbse_logger.setPreface(ctp::logERROR, (format("\nGWBSE ERR ...")).str());
          gwbse_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...")).str());
          gwbse_logger.setPreface(ctp::logDEBUG, (format("\nGWBSE DBG ...")).str());
          GWBSE gwbse = GWBSE(orbitalsAB);
          gwbse.setLogger(&gwbse_logger);
          gwbse.Initialize(_gwbse_options);
          gwbse.Evaluate();
           WriteLoggerToFile(work_dir + "/gwbse.log", gwbse_logger);
        } catch (std::runtime_error& error) {
          std::string errormessage(error.what());
          SetJobToFailed(jres, pLog, errormessage);
          return jres;
        }

      } // end of excited state calculation, exciton data is in _orbitalsAB

      // calculate the coupling


      if (_do_bsecoupling) {
        CTP_LOG(ctp::logDEBUG, *pLog) << "Running BSECoupling" << flush;
        BSECoupling bsecoupling;
        // orbitals must be loaded from a file
        if (!_do_gwbse) {
          try {
            orbitalsAB.ReadFromCpt(orbFileAB);
          } catch (std::runtime_error& error) {
            SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileAB);
            return jres;
          }
        }

        Orbitals orbitalsB;
        Orbitals orbitalsA;

        try {
          orbitalsA.ReadFromCpt(orbFileA);
        } catch (std::runtime_error& error) {
          SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileA);
          return jres;
        }

        try {
          orbitalsB.ReadFromCpt(orbFileB);
        } catch (std::runtime_error& error) {
          SetJobToFailed(jres, pLog, "Do input: failed loading orbitals from " + orbFileB);
          return jres;
        }

        try {
          ctp::Logger bsecoupling_logger(ctp::logDEBUG);
          bsecoupling_logger.setMultithreading(false);
          bsecoupling_logger.setPreface(ctp::logINFO, (format("\nGWBSE INF ...")).str());
          bsecoupling_logger.setPreface(ctp::logERROR, (format("\nGWBSE ERR ...")).str());
          bsecoupling_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...")).str());
          bsecoupling_logger.setPreface(ctp::logDEBUG, (format("\nGWBSE DBG ...")).str());
          bsecoupling.setLogger(&bsecoupling_logger);
          bsecoupling.Initialize(_bsecoupling_options);
          bsecoupling.CalculateCouplings(orbitalsA, orbitalsB, orbitalsAB);
          bsecoupling.Addoutput(job_output, orbitalsA, orbitalsB);
          WriteLoggerToFile(work_dir + "/bsecoupling.log", bsecoupling_logger);
        } catch (std::runtime_error& error) {
          std::string errormessage(error.what());
          SetJobToFailed(jres, pLog, errormessage);
          return jres;
        }

      }

      tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");
      stringstream sout;
      sout << iomXML << _job_summary;
      CTP_LOG(ctp::logINFO, *pLog) << ctp::TimeStamp() << " Finished evaluating pair " << ID_A << ":" << ID_B << flush;
      if (_store_dft || _store_singlets || _store_triplets || _store_ehint) {
        boost::filesystem::create_directories(orb_dir);
        CTP_LOG(ctp::logDEBUG, *pLog) << "Saving orbitals to " << orbFileAB << flush;
        if (!_store_dft) {
          orbitalsAB.AOVxc().resize(0, 0);
          orbitalsAB.MOCoefficients().resize(0, 0);
        }
        if (!_store_singlets) {
          orbitalsAB.BSESingletCoefficients().resize(0, 0);
          orbitalsAB.BSESingletEnergies().resize(0, 0);
        }
        if (!_store_triplets) {
          orbitalsAB.BSETripletCoefficients().resize(0, 0);
          orbitalsAB.BSETripletEnergies().resize(0, 0);
        }
        if (!_store_ehint) {
          orbitalsAB.eh_t().resize(0, 0);
          orbitalsAB.eh_s().resize(0, 0);
        }
        orbitalsAB.WriteToCpt(orbFileAB);
      } else {
        CTP_LOG(ctp::logDEBUG, *pLog) << "Orb file is not saved according to options " << flush;
      }

      jres.setOutput(_job_summary);
      jres.setStatus(ctp::Job::COMPLETE);

      return jres;
    }

    void IQM::WriteJobFile(ctp::Topology * top) {

      cout << endl << "... ... Writing job file " << flush;
      std::ofstream ofs;
      ofs.open(_jobfile.c_str(), std::ofstream::out);
      if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);

      ctp::QMNBList::iterator pit;
      ctp::QMNBList &nblist = top->NBList();

      int jobCount = 0;
      if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
      }

      ofs << "<jobs>" << endl;
      string tag = "";

      for (ctp::QMPair* pair : nblist) {
        if (pair->getType() == ctp::QMPair::Excitoncl) {
          continue;
        }
        int id1 = pair->Seg1()->getId();
        string name1 = pair->Seg1()->getName();
        int id2 = pair->Seg2()->getId();
        string name2 = pair->Seg2()->getName();
        int id = pair->getId();
        Property Input;
        Property &pInput = Input.add("input", "");
        Property &pSegmentA = pInput.add("segment", boost::lexical_cast<string>(id1));
        pSegmentA.setAttribute<string>("type", name1);
        pSegmentA.setAttribute<int>("id", id1);
        Property & pSegmentB = pInput.add("segment", boost::lexical_cast<string>(id2));
        pSegmentB.setAttribute<string>("type", name2);
        pSegmentB.setAttribute<int>("id", id2);
        ctp::Job job(id, tag, Input, ctp::Job::AVAILABLE);
        job.ToStream(ofs, "xml");
      }
      // CLOSE STREAM
      ofs << "</jobs>" << endl;
      ofs.close();
      cout << endl << "... ... In total " << jobCount << " jobs" << flush;
      return;
    }

    double IQM::GetDFTCouplingFromProp(tools::Property& dftprop, int stateA, int stateB) {
      double J = 0;
      double found=false;
      for (Property* state : dftprop.Select("coupling")) {
        int state1 = state->getAttribute<int>("levelA");
        int state2 = state->getAttribute<int>("levelB");
        if (state1 == stateA && state2 == stateB) {
          J = state->getAttribute<double>("j");
          found=true;
          break;
        }

      }
     if(found){
      return J*J;
      }else{
        return -1;
      }
    }

    double IQM::GetBSECouplingFromProp(tools::Property& bseprop,const QMState& stateA,const QMState& stateB) {
      double J = 0;
      std::string algorithm=bseprop.getAttribute<std::string>("algorithm");
      double found=false;
      for (Property* state : bseprop.Select("coupling")) {
        QMState state1;
        state1.FromString(state->getAttribute<std::string>("stateA"));
        QMState state2;
        state2.FromString(state->getAttribute<std::string>("stateB"));
        if (state1 == stateA && state2 == stateB) {
          J = state->getAttribute<double>(algorithm);
          found=true;
          break;
        }
      }
      if(found){
      return J*J;
      }else{
        return -1;
      }
    }
    
    
    QMState IQM::GetElementFromMap(const std::map<std::string, QMState>& elementmap,const std::string& elementname )const{
      QMState state;
      try{
        state = elementmap.at(elementname);
      }
      catch (std::out_of_range& error) {
        std::string errormessage="Map does not have segment of type: "+elementname;
        errormessage+="\n segments in map are:";
        for(const auto& s:elementmap){
         errormessage+="\n\t"+s.first;
        }
        throw std::runtime_error(errormessage);
      }
      return state;
    }

    void IQM::ReadJobFile(ctp::Topology * top) {
      // gets the neighborlist from the topology
      ctp::QMNBList &nblist = top->NBList();
      int number_of_pairs = nblist.size();
      int dft_h = 0;
      int dft_e = 0;
      int bse_s = 0;
      int bse_t = 0;
      int incomplete_jobs = 0;
      ctp::Logger log;
      log.setReportLevel(ctp::logINFO);
        
      Property xml;
      // load the QC results in a vector indexed by the pair ID
      load_property_from_xml(xml, _jobfile);
      list<Property*> jobProps = xml.Select("jobs.job");
      vector<Property*> records = std::vector<Property*>(nblist.size() + 1, NULL);

      // loop over all jobs = pair records in the job file
      for (Property* job : jobProps) {
        if (job->exists("status")) {
          if (job->get("status").as<std::string>() != "COMPLETE" || !job->exists("output")) {
            incomplete_jobs++;
            continue;
          }
        }

        // get the output records
        Property poutput = job->get("output");
        // job file is stupid, because segment ids are only in input have to get them out l
        list<Property*> segmentprobs = job->Select("input.segment");
        vector<int> id;
        for (Property* segment : segmentprobs) {
          id.push_back(segment->getAttribute<int>("id"));
        }
        if (id.size() != 2) throw std::runtime_error("Getting pair ids from jobfile failed, check jobfile.");

        double idA = id[0];
        double idB = id[1];

        // segments which correspond to these ids           
        ctp::Segment *segA = top->getSegment(idA);
        ctp::Segment *segB = top->getSegment(idB);
        // pair that corresponds to the two segments
        ctp::QMPair *qmp = nblist.FindPair(segA, segB);
        // output using logger
        
        if (qmp == NULL) { // there is no pair in the neighbor list with this name
          CTP_LOG_SAVE(ctp::logINFO, log) << "No pair " << idA << ":" << idB << " found in the neighbor list. Ignoring" << flush;
        } else {
          records[qmp->getId()] = &(job->get("output"));
        }

      } // finished loading from the file

      for (ctp::QMPair *pair:top->NBList()) {

        if (records[ pair->getId() ] == NULL) continue; //skip pairs which are not in the jobfile

        ctp::Segment* segmentA = pair->Seg1();
        ctp::Segment* segmentB = pair->Seg2();

        ctp::QMPair::PairType ptype = pair->getType();
        if (ptype != ctp::QMPair::PairType::Hopping
                && ptype != ctp::QMPair::PairType::SuperExchangeAndHopping) {
          cout << "WARNING Pair " << pair->getId() << " is not of any of the "
                  "Hopping or SuperExchangeAndHopping type. Skipping pair" << flush;
          continue;
        }
        
        Property* pair_property = records[ pair->getId() ];

        if (pair_property->exists("dftcoupling")) {
          tools::Property& dftprop = pair_property->get("dftcoupling");
          int homoA = dftprop.getAttribute<int>("homoA");
          int homoB = dftprop.getAttribute<int>("homoB");
          QMStateType hole=QMStateType(QMStateType::Hole);
          if (dftprop.exists(hole.ToLongString())) {
            tools::Property& holes = dftprop.get(hole.ToLongString());
            QMState stateA = GetElementFromMap(_hole_levels,segmentA->getName());
            QMState stateB = GetElementFromMap(_hole_levels,segmentB->getName());
            int levelA = homoA - stateA.Index(); //h1 is is homo;
            int levelB = homoB - stateB.Index();
            double J2 = GetDFTCouplingFromProp(holes, levelA, levelB);
            if(J2>0){
              pair->setJeff2(J2, 1);
              pair->setIsPathCarrier(true, 1);
              dft_h++;
            }
          }
          QMStateType electron=QMStateType(QMStateType::Electron);
          if (dftprop.exists(electron.ToLongString())) {
            tools::Property& electrons = dftprop.get(electron.ToLongString());
            QMState stateA = GetElementFromMap(_electron_levels,segmentA->getName());
            QMState stateB = GetElementFromMap(_electron_levels,segmentB->getName());
            int levelA = homoA +1+stateA.Index(); //e1 is lumo;
            int levelB = homoB +1+stateB.Index();
            double J2 = GetDFTCouplingFromProp(electrons, levelA, levelB);
            if(J2>0){
              pair->setJeff2(J2, -1);
              pair->setIsPathCarrier(true, -1);
              dft_e++;
            }
          }
        }
        if (pair_property->exists("bsecoupling")) {
          tools::Property& bseprop = pair_property->get("bsecoupling");
          QMStateType singlet=QMStateType(QMStateType::Singlet);
          if (bseprop.exists(singlet.ToLongString())) {
            tools::Property& singlets = bseprop.get(singlet.ToLongString());
            QMState stateA = GetElementFromMap(_singlet_levels,segmentA->getName());
            QMState stateB = GetElementFromMap(_singlet_levels,segmentB->getName());
            double J2 = GetBSECouplingFromProp(singlets, stateA, stateB);
            if(J2>0){
              pair->setJeff2(J2, 2);
              pair->setIsPathCarrier(true, 2);
              bse_s++;
            }
          }
          QMStateType triplet=QMStateType(QMStateType::Triplet);
          if (bseprop.exists(triplet.ToLongString())) {
            tools::Property& triplets = bseprop.get(triplet.ToLongString());
            QMState stateA = GetElementFromMap(_triplet_levels,segmentA->getName());
            QMState stateB = GetElementFromMap(_triplet_levels,segmentB->getName());
            double J2 = GetBSECouplingFromProp(triplets, stateA, stateB);
            if(J2>0){
              pair->setJeff2(J2, 3);
              pair->setIsPathCarrier(true, 3);
              bse_t++;
            }
          }
        }

      }

      CTP_LOG_SAVE(ctp::logINFO, log) << "Pairs [total:updated(e,h,s,t)] "
              << number_of_pairs << ":(" << dft_e << ","<< dft_h << "," << bse_s<< "," << bse_t
              << ") Incomplete jobs: " << incomplete_jobs << flush;
      cout<<std::endl;
      cout << log;
      return;
    }

  }
};
