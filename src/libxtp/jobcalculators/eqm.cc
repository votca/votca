/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#include "votca/xtp/segmentmapper.h"
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/esp2multipole.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
namespace xtp {

void EQM::Initialize(tools::Property& options) {

  _do_dft_input = false;
  _do_dft_run = false;
  _do_dft_parse = false;
  _do_gwbse = false;
  _do_esp = false;

  ParseOptionsXML(options);
  QMPackageFactory::RegisterAll();
}

void EQM::ParseOptionsXML(tools::Property& options) {

  _maverick = (_nThreads == 1) ? true : false;
  std::string key = "options." + Identify();
  // job tasks
  std::string _tasks_string = options.get(key + ".tasks").as<std::string>();
  if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
  if (_tasks_string.find("dft") != std::string::npos) _do_dft_run = true;
  if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
  if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;
  if (_tasks_string.find("esp") != std::string::npos) _do_esp = true;

  // options for gwbse
  key = "options." + Identify();
  std::string _gwbse_xml =
      options.get(key + ".gwbse_options").as<std::string>();
  load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());

  // options for dft package
  std::string _package_xml = options.get(key + ".dftpackage").as<std::string>();
  load_property_from_xml(_package_options, _package_xml.c_str());
  key = "package";
  _package = _package_options.get(key + ".name").as<std::string>();

  // options for esp/partialcharges
  if (_do_esp) {
    key = "options." + Identify();
    std::string _esp_xml = options.get(key + ".esp_options").as<std::string>();
    load_property_from_xml(_esp_options, _esp_xml.c_str());
  }
  _jobfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".job_file");
  _mapfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".map_file");
}

void EQM::WriteJobFile(Topology& top) {

  std::cout << std::endl << "... ... Writing job file: " << std::flush;
  std::ofstream ofs;
  ofs.open(_jobfile, std::ofstream::out);
  if (!ofs.is_open())
    throw std::runtime_error("\nERROR: bad file handle: " + _jobfile);
  ofs << "<jobs>" << std::endl;
  int jobCount = 0;

  const std::vector<Segment>& segments = top.Segments();
  for (const Segment& segment : segments) {
    int id = ++jobCount;
    std::string tag = "";
    tools::Property Input;
    tools::Property& pInput = Input.add("input", "");
    tools::Property& pSegment =
        pInput.add("segment", (format("%1$s") % segment.getId()).str());
    pSegment.setAttribute<std::string>("type", segment.getName());
    pSegment.setAttribute<int>("id", segment.getId());
    Job job(id, tag, Input, Job::AVAILABLE);
    job.ToStream(ofs, "xml");
  }

  ofs << "</jobs>" << std::endl;
  ofs.close();

  std::cout << jobCount << " jobs" << std::flush;
}

void EQM::SetJobToFailed(Job::JobResult& jres, Logger& pLog,
                         const std::string& errormessage) {
  XTP_LOG_SAVE(logERROR, pLog) << errormessage << std::flush;
  std::cout << pLog;
  jres.setError(errormessage);
  jres.setStatus(Job::FAILED);
}

void EQM::WriteLoggerToFile(const std::string& logfile, Logger& logger) {
  std::ofstream ofs;
  ofs.open(logfile, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + logfile);
  }
  ofs << logger << std::endl;
  ofs.close();
}
Job::JobResult EQM::EvalJob(Topology& top, Job& job, QMThread& opThread) {

  Orbitals orbitals;
  Job::JobResult jres = Job::JobResult();
  tools::Property _job_input = job.getInput();
  std::vector<tools::Property*> lSegments = _job_input.Select("segment");
  int segId = lSegments.front()->getAttribute<int>("id");
  std::string segType = lSegments.front()->getAttribute<std::string>("type");
  std::string qmgeo_state = "n";
  if (lSegments.front()->exists("qm_geometry")) {
    qmgeo_state = lSegments.front()->getAttribute<std::string>("qm_geometry");
  }

  QMState state(qmgeo_state);
  Segment& seg = top.getSegment(segId);

  Logger& pLog = opThread.getLogger();
  QMMapper mapper(pLog);
  mapper.LoadMappingFile(_mapfile);
  orbitals.QMAtoms() = mapper.map(seg, state);
  XTP_LOG_SAVE(logINFO, pLog)
      << TimeStamp() << " Evaluating site " << seg.getId() << std::flush;

  // directories and files
  boost::filesystem::path arg_path;
  std::string eqm_work_dir = "OR_FILES";
  std::string frame_dir =
      "frame_" + boost::lexical_cast<std::string>(top.getStep());
  std::string orb_file =
      (format("%1%_%2%%3%") % "molecule" % segId % ".orb").str();
  std::string mol_dir = (format("%1%%2%%3%") % "molecule" % "_" % segId).str();
  std::string package_append = _package + "_" + Identify();
  std::string work_dir =
      (arg_path / eqm_work_dir / package_append / frame_dir / mol_dir).c_str();

  tools::Property job_summary;
  tools::Property& output_summary = job_summary.add("output", "");
  tools::Property& segment_summary = output_summary.add("segment", "");
  std::string segName = seg.getName();
  segId = seg.getId();
  segment_summary.setAttribute("id", segId);
  segment_summary.setAttribute("type", segName);
  if (_do_dft_input || _do_dft_run || _do_dft_parse) {
    XTP_LOG_SAVE(logDEBUG, pLog) << "Running DFT" << std::flush;
    Logger dft_logger(logDEBUG);
    dft_logger.setMultithreading(false);
    dft_logger.setPreface(logINFO, (format("\nDFT INF ...")).str());
    dft_logger.setPreface(logERROR, (format("\nDFT ERR ...")).str());
    dft_logger.setPreface(logWARNING, (format("\nDFT WAR ...")).str());
    dft_logger.setPreface(logDEBUG, (format("\nDFT DBG ...")).str());

    std::unique_ptr<QMPackage> qmpackage =
        std::unique_ptr<QMPackage>(QMPackages().Create(_package));
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
      run_dft_status = qmpackage->Run();
      if (!run_dft_status) {
        std::string output = "DFT run failed";
        SetJobToFailed(jres, pLog, output);
        return jres;
      }
    }

    // parse the log/orbitals files
    bool parse_log_status = false;
    bool parse_orbitals_status = false;
    if (_do_dft_parse) {
      parse_log_status = qmpackage->ParseLogFile(orbitals);
      if (!parse_log_status) {
        std::string output = "log incomplete; ";
        SetJobToFailed(jres, pLog, output);
        return jres;
      }
      parse_orbitals_status = qmpackage->ParseMOsFile(orbitals);
      if (!parse_orbitals_status) {
        std::string output = "orbfile failed; ";
        SetJobToFailed(jres, pLog, output);
        return jres;
      }
    }  // end of the parse orbitals/log
    qmpackage->CleanUp();
    WriteLoggerToFile(work_dir + "/dft.log", dft_logger);
  }

  if (!_do_dft_parse) {
    // load the DFT data from serialized orbitals object
    std::string ORB_FILE =
        eqm_work_dir + "/molecules/" + frame_dir + "/" + orb_file;
    XTP_LOG_SAVE(logDEBUG, pLog)
        << TimeStamp() << " Loading DFT data from " << ORB_FILE << std::flush;
    orbitals.ReadFromCpt(ORB_FILE);
  }

  if (_do_gwbse) {
    XTP_LOG_SAVE(logDEBUG, pLog) << "Running GWBSE" << std::flush;
    try {
      GWBSE gwbse = GWBSE(orbitals);
      Logger gwbse_logger(logDEBUG);
      gwbse_logger.setMultithreading(false);
      gwbse_logger.setPreface(logINFO, (format("\nGWBSE INF ...")).str());
      gwbse_logger.setPreface(logERROR, (format("\nGWBSE ERR ...")).str());
      gwbse_logger.setPreface(logWARNING, (format("\nGWBSE WAR ...")).str());
      gwbse_logger.setPreface(logDEBUG, (format("\nGWBSE DBG ...")).str());
      gwbse.setLogger(&gwbse_logger);
      gwbse.Initialize(_gwbse_options);
      gwbse.Evaluate();
      gwbse.addoutput(segment_summary);
      WriteLoggerToFile(work_dir + "/gwbse.log", gwbse_logger);
    } catch (std::runtime_error& error) {
      std::string errormessage(error.what());
      SetJobToFailed(jres, pLog, "GWBSE:" + errormessage);
      return jres;
    }
  }

  if (_do_esp) {
    XTP_LOG_SAVE(logDEBUG, pLog) << "Running ESPFIT" << std::flush;
    try {
      Esp2multipole esp2multipole = Esp2multipole(pLog);
      esp2multipole.Initialize(_esp_options);
      std::string ESPDIR =
          "MP_FILES/" + frame_dir + "/" + esp2multipole.GetStateString();
      esp2multipole.Extractingcharges(orbitals);
      std::string mps_file = (format("%1%_%2%_%3%.mps") % segType % segId %
                  esp2multipole.GetStateString())
                     .str();
      boost::filesystem::create_directories(ESPDIR);
      esp2multipole.WritetoFile(ESPDIR + "/" + mps_file, orbitals);
      XTP_LOG_SAVE(logDEBUG, pLog)
          << "Written charges to " << (ESPDIR + "/" + mps_file).c_str()
          << std::flush;
      segment_summary.add("partialcharges", (ESPDIR + "/" + mps_file).c_str());
    } catch (std::runtime_error& error) {
      std::string errormessage(error.what());
      SetJobToFailed(jres, pLog, "ESPFIT:" + errormessage);
      return jres;
    }
  }
  XTP_LOG_SAVE(logINFO, pLog) << TimeStamp() << " Finished evaluating site "
                              << seg.getId() << std::flush;

  if (_do_dft_parse || _do_gwbse) {
    XTP_LOG_SAVE(logDEBUG, pLog) << "Saving data to " << orb_file << std::flush;
    std::string DIR = eqm_work_dir + "/molecules/" + frame_dir;
    boost::filesystem::create_directories(DIR);
    std::string ORBFILE = DIR + "/" + orb_file;
    orbitals.WriteToCpt(ORBFILE);
  }

  // output of the JOB
  jres.setOutput(job_summary);
  jres.setStatus(Job::COMPLETE);

  return jres;
}
}  // namespace xtp
};  // namespace votca
