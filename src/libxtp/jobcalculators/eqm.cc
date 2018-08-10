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


#include "eqm.h"
#include <votca/xtp/esp2multipole.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>



#include <boost/math/constants/constants.hpp>


using boost::format;
using namespace boost::filesystem;

namespace votca {
  namespace xtp {

    void EQM::Initialize(Property *options) {

      _do_dft_input = false;
      _do_dft_run = false;
      _do_dft_parse = false;
      _do_gwbse = false;
      _do_esp = false;

      ParseOptionsXML(options);
      QMPackageFactory::RegisterAll();
    }

    void EQM::ParseOptionsXML(Property* options) {

      _maverick = (_nThreads == 1) ? true : false;
      string key = "options." + Identify();
      // job tasks
      string _tasks_string = options->get(key + ".tasks").as<string> ();
      if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
      if (_tasks_string.find("dft") != std::string::npos) _do_dft_run = true;
      if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
      if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;
      if (_tasks_string.find("esp") != std::string::npos) _do_esp = true;

      key = "options." + Identify();

      if (options->exists(key + ".job_file")) {
        _jobfile = options->get(key + ".job_file").as<string>();
      } else {
        throw std::runtime_error("Job-file not set. Abort.");
      }

      // options for gwbse
      key = "options." + Identify();
      string _gwbse_xml = options->get(key + ".gwbse_options").as<string> ();
      load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());

      // options for dft package
      string _package_xml = options->get(key + ".dftpackage").as<string> ();
      load_property_from_xml(_package_options, _package_xml.c_str());
      key = "package";
      _package = _package_options.get(key + ".name").as<string> ();

      //options for esp/partialcharges
      if (_do_esp) {
        key = "options." + Identify();
        string _esp_xml = options->get(key + ".esp_options").as<string> ();
        load_property_from_xml(_esp_options, _esp_xml.c_str());
      }

    }

    void EQM::WriteJobFile(ctp::Topology *top) {

      cout << endl << "... ... Writing job file: " << flush;
      std::ofstream ofs;
      ofs.open(_jobfile.c_str(), std::ofstream::out);
      if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);
      ofs << "<jobs>" << endl;
      int jobCount = 0;

      std::vector<ctp::Segment*> segments = top->Segments();
      for (ctp::Segment* segment : segments) {
        int id = ++jobCount;
        string tag = "";
        Property Input;
        Property &pInput = Input.add("input", "");
        Property &pSegment = pInput.add("segment", (format("%1$s") % segment->getId()).str());
        pSegment.setAttribute<string>("type", segment->getName());
        pSegment.setAttribute<int>("id", segment->getId());
        ctp::Job job(id, tag, Input, ctp::Job::AVAILABLE);
        job.ToStream(ofs, "xml");
      }

      ofs << "</jobs>" << endl;
      ofs.close();

      cout << jobCount << " jobs" << flush;

    }

    void EQM::WriteLoggerToFile(const string& logfile, ctp::Logger& logger){
      std::ofstream ofs;
      ofs.open(logfile.c_str(), std::ofstream::out);
      if (!ofs.is_open()) {
        throw runtime_error("Bad file handle: " + logfile);
      }
      ofs << logger << endl;
      ofs.close();
    }

    ctp::Job::JobResult EQM::EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *opThread) {

      Orbitals orbitals;
      ctp::Job::JobResult jres = ctp::Job::JobResult();
      Property _job_input = job->getInput();
      list<Property*> lSegments = _job_input.Select("segment");
      vector < ctp::Segment* > segments;
      int segId = lSegments.front()->getAttribute<int>("id");
      string segType = lSegments.front()->getAttribute<string>("type");
      ctp::Segment *seg = top->getSegment(segId);
      assert(seg->getName() == segType);
      segments.push_back(seg);
      QMInterface interface;
      orbitals.QMAtoms() = interface.Convert(segments);

      ctp::Logger* pLog = opThread->getLogger();

      CTP_LOG(ctp::logINFO, *pLog) << ctp::TimeStamp() << " Evaluating site " << seg->getId() << flush;

      // directories and files
      path arg_path;
      string eqm_work_dir = "OR_FILES";
      string frame_dir = "frame_" + boost::lexical_cast<string>(top->getDatabaseId());
      string orb_file = (format("%1%_%2%%3%") % "molecule" % segId % ".orb").str();
      string mol_dir = (format("%1%%2%%3%") % "molecule" % "_" % segId).str();
      string package_append = _package + "_"+Identify();
      string work_dir = (arg_path / eqm_work_dir / package_append / frame_dir / mol_dir).c_str();

      Property job_summary;
      Property &output_summary = job_summary.add("output", "");
      Property &segment_summary = output_summary.add("segment", "");
      string segName = seg->getName();
      segId = seg->getId();
      segment_summary.setAttribute("id", segId);
      segment_summary.setAttribute("type", segName);
      if (_do_dft_input || _do_dft_run || _do_dft_parse) {
        
        ctp::Logger dft_logger(ctp::logDEBUG);
        dft_logger.setMultithreading(false);
        dft_logger.setPreface(ctp::logINFO, (format("\nDFT INF ...")).str());
        dft_logger.setPreface(ctp::logERROR, (format("\nDFT ERR ...")).str());
        dft_logger.setPreface(ctp::logWARNING, (format("\nDFT WAR ...")).str());
        dft_logger.setPreface(ctp::logDEBUG, (format("\nDFT DBG ...")).str());

        QMPackage *qmpackage = QMPackages().Create(_package);
        qmpackage->setLog(&dft_logger);
        qmpackage->setRunDir(work_dir);
        qmpackage->Initialize(_package_options);

        // create input for DFT
        if (_do_dft_input) {
          boost::filesystem::create_directories(work_dir);
          qmpackage->WriteInputFile(orbitals);
        }

        bool run_dft_status = false;
        if (_do_dft_run) {
          run_dft_status = qmpackage->Run(orbitals);
          if (!run_dft_status) {
            string output = "run failed; ";
            CTP_LOG(ctp::logERROR, *pLog) << qmpackage->getPackageName() << " run failed" << flush;
            cout << *pLog;
            jres.setOutput(output);
            jres.setStatus(ctp::Job::FAILED);
            delete qmpackage;
            return jres;
          }
        }

        // parse the log/orbitals files
        bool parse_log_status = false;
        bool parse_orbitals_status = false;
        if (_do_dft_parse) {
          parse_log_status = qmpackage->ParseLogFile(orbitals);
          if (!parse_log_status) {
            string output = "log incomplete; ";
            CTP_LOG(ctp::logERROR, *pLog) << "LOG parsing failed" << flush;
            cout << *pLog;
            jres.setOutput(output);
            jres.setStatus(ctp::Job::FAILED);
            delete qmpackage;
            return jres;
          }
          parse_orbitals_status = qmpackage->ParseOrbitalsFile(orbitals);
          if (!parse_orbitals_status) {
            string output = "orbfile failed; ";
            CTP_LOG(ctp::logERROR, *pLog) << "Orbitals parsing failed" << flush;
            cout << *pLog;
            jres.setOutput(output);
            jres.setStatus(ctp::Job::FAILED);
            delete qmpackage;
            return jres;
          }
        }// end of the parse orbitals/log
        qmpackage->CleanUp();
        delete qmpackage;
         WriteLoggerToFile(work_dir + "/dft.log",dft_logger);
      }

      if (!_do_dft_parse) {
        // load the DFT data from serialized orbitals object
        string ORB_FILE =eqm_work_dir + "/molecules/" + frame_dir+ "/" + orb_file;
        CTP_LOG(ctp::logDEBUG, *pLog) << ctp::TimeStamp() << " Loading DFT data from " << ORB_FILE << flush;
        orbitals.ReadFromCpt(ORB_FILE);
      }

      if (_do_gwbse) {

        GWBSE gwbse = GWBSE(orbitals);
        
        // define own logger for GW-BSE that is written into a runFolder logfile
        ctp::Logger gwbse_logger(ctp::logDEBUG);
        gwbse_logger.setMultithreading(false);
        gwbse_logger.setPreface(ctp::logINFO, (format("\nGWBSE INF ...")).str());
        gwbse_logger.setPreface(ctp::logERROR, (format("\nGWBSE ERR ...")).str());
        gwbse_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...")).str());
        gwbse_logger.setPreface(ctp::logDEBUG, (format("\nGWBSE DBG ...")).str());
        gwbse.setLogger(&gwbse_logger);
        gwbse.Initialize(_gwbse_options);
        gwbse.Evaluate();
        gwbse.addoutput(segment_summary);
        // write logger to log file
        WriteLoggerToFile(work_dir + "/gwbse.log", gwbse_logger);

      }

      if (_do_esp) {
        string mps_file = "";
        Esp2multipole esp2multipole = Esp2multipole(pLog);
        esp2multipole.Initialize(_esp_options);
        string ESPDIR = "MP_FILES/" + frame_dir + "/" + esp2multipole.GetIdentifier();
        esp2multipole.Extractingcharges(orbitals);
        mps_file = (format("%1%_%2%_%3%.mps") % segType % segId % esp2multipole.GetIdentifier()).str();
        boost::filesystem::create_directories(ESPDIR);
        esp2multipole.WritetoFile((ESPDIR + "/" + mps_file).c_str(), Identify());
        CTP_LOG(ctp::logDEBUG, *pLog) << "Written charges to " << (ESPDIR + "/" + mps_file).c_str() << flush;
        segment_summary.add("partialcharges", (ESPDIR + "/" + mps_file).c_str());
      }
      CTP_LOG(ctp::logINFO, *pLog) << ctp::TimeStamp() << " Finished evaluating site " << seg->getId() << flush;

      if (_do_dft_parse || _do_gwbse) {
        CTP_LOG(ctp::logDEBUG, *pLog) << "Saving data to " << orb_file << flush;
        string DIR = eqm_work_dir + "/molecules/" + frame_dir;
        boost::filesystem::create_directories(DIR);
        string ORBFILE = DIR + "/" + orb_file;
        orbitals.WriteToCpt(ORBFILE);
      }

      // output of the JOB 
      jres.setOutput(job_summary);
      jres.setStatus(ctp::Job::COMPLETE);

      // dump the LOG

      return jres;
    }

  }


};
