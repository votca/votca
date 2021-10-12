/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Third party includes
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

// Local VOTCA includes
#include "votca/xtp/esp2multipole.h"
#include "votca/xtp/segmentmapper.h"

// Local private VOTCA includes
#include "eqm.h"

using boost::format;
using namespace boost::filesystem;

namespace votca {
namespace xtp {

void EQM::ParseSpecificOptions(const tools::Property& options) {

  QMPackageFactory::RegisterAll();

  // job tasks
  std::string tasks_string_ = options.get(".tasks").as<std::string>();

  if (tasks_string_.find("input") != std::string::npos) {
    do_dft_input_ = true;
  }
  if (tasks_string_.find("dft") != std::string::npos) {
    do_dft_run_ = true;
  }
  if (tasks_string_.find("parse") != std::string::npos) {
    do_dft_parse_ = true;
  }
  if (tasks_string_.find("gwbse") != std::string::npos) {
    do_gwbse_ = true;
  }
  if (tasks_string_.find("esp") != std::string::npos) {
    do_esp_ = true;
  }

  // set the basis sets and functional for DFT and GWBSE
  gwbse_options_ = options.get("gwbse");
  package_options_ = options.get("dftpackage");
  esp_options_ = options.get(".esp_options");
}

void EQM::WriteJobFile(const Topology& top) {

  std::cout << "\n... ... Writing job file: " << jobfile_ << std::flush;
  std::ofstream ofs;
  ofs.open(jobfile_, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("\nERROR: bad file handle: " + jobfile_);
  }
  ofs << "<jobs>" << std::endl;
  Index jobCount = 0;

  const std::vector<Segment>& segments = top.Segments();
  for (const Segment& segment : segments) {
    Index id = segment.getId();
    std::string tag = "";
    tools::Property Input;
    tools::Property& pInput = Input.add("input", "");
    tools::Property& pSegment =
        pInput.add("segment", (format("%1$s") % segment.getId()).str());
    pSegment.setAttribute<std::string>("type", segment.getType());
    pSegment.setAttribute<Index>("id", segment.getId());
    Job job(id, tag, Input, Job::AVAILABLE);
    job.ToStream(ofs);
    jobCount++;
  }

  ofs << "</jobs>" << std::endl;
  ofs.close();

  std::cout << " with " << jobCount << " jobs" << std::flush;
}

void EQM::SetJobToFailed(Job::JobResult& jres, Logger& pLog,
                         const std::string& errormessage) {
  XTP_LOG(Log::error, pLog) << errormessage << std::flush;
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

Job::JobResult EQM::EvalJob(const Topology& top, Job& job, QMThread& opThread) {
  Orbitals orbitals;
  Job::JobResult jres = Job::JobResult();
  tools::Property job_input_ = job.getInput();
  std::vector<tools::Property*> lSegments = job_input_.Select("segment");
  Index segId = lSegments.front()->getAttribute<Index>("id");
  std::string segType = lSegments.front()->getAttribute<std::string>("type");
  std::string qmgeo_state = "n";
  if (lSegments.front()->exists("qm_geometry")) {
    qmgeo_state = lSegments.front()->getAttribute<std::string>("qm_geometry");
  }

  QMState state(qmgeo_state);
  const Segment& seg = top.getSegment(segId);

  Logger& pLog = opThread.getLogger();
  QMMapper mapper(pLog);
  mapper.LoadMappingFile(mapfile_);
  orbitals.QMAtoms() = mapper.map(seg, state);
  XTP_LOG(Log::error, pLog)
      << TimeStamp() << " Evaluating site " << seg.getId() << std::flush;

  // directories and files
  boost::filesystem::path arg_path;
  std::string eqm_work_dir = "OR_FILES";
  std::string frame_dir =
      "frame_" + boost::lexical_cast<std::string>(top.getStep());
  std::string orb_file =
      (format("%1%_%2%%3%") % "molecule" % segId % ".orb").str();
  std::string mol_dir = (format("%1%%2%%3%") % "molecule" % "_" % segId).str();
  std::string package_append = "workdir_" + Identify();
  std::string work_dir =
      (arg_path / eqm_work_dir / package_append / frame_dir / mol_dir)
          .generic_string();

  tools::Property job_summary;
  tools::Property& output_summary = job_summary.add("output", "");
  tools::Property& segment_summary = output_summary.add("segment", "");
  std::string segName = seg.getType();
  segId = seg.getId();
  segment_summary.setAttribute("id", segId);
  segment_summary.setAttribute("type", segName);
  if (do_dft_input_ || do_dft_run_ || do_dft_parse_) {
    XTP_LOG(Log::error, pLog) << "Running DFT" << std::flush;
    Logger dft_logger(votca::Log::current_level);
    dft_logger.setMultithreading(false);
    dft_logger.setPreface(Log::info, (format("\nDFT INF ...")).str());
    dft_logger.setPreface(Log::error, (format("\nDFT ERR ...")).str());
    dft_logger.setPreface(Log::warning, (format("\nDFT WAR ...")).str());
    dft_logger.setPreface(Log::debug, (format("\nDFT DBG ...")).str());
    std::string package = package_options_.get(".name").as<std::string>();
    std::unique_ptr<QMPackage> qmpackage =
        QMPackageFactory::QMPackages().Create(package);
    qmpackage->setLog(&dft_logger);
    qmpackage->setRunDir(work_dir);
    qmpackage->Initialize(package_options_);

    // create input for DFT
    if (do_dft_input_) {
      boost::filesystem::create_directories(work_dir);
      qmpackage->WriteInputFile(orbitals);
    }

    if (do_dft_run_) {
      bool run_dft_status = qmpackage->Run(orbitals);
      if (!run_dft_status) {
        std::string output = "DFT run failed";
        SetJobToFailed(jres, pLog, output);
        return jres;
      }
      // additionally copy *.gbw files for orca (-> iqm guess)
      if (qmpackage->getPackageName() == "orca") {
        std::string DIR = eqm_work_dir + "/molecules/" + frame_dir;
        boost::filesystem::create_directories(DIR);
        std::string gbw_file =
            (format("%1%_%2%%3%") % "molecule" % segId % ".gbw").str();
        std::string GBWFILE = DIR + "/" + gbw_file;
        XTP_LOG(Log::error, pLog)
            << "Copying MO data to " << gbw_file << std::flush;
        std::string GBWFILE_workdir = work_dir + "/system.gbw";
        boost::filesystem::copy_file(
            GBWFILE_workdir, GBWFILE,
            boost::filesystem::copy_option::overwrite_if_exists);
      }
    }
    WriteLoggerToFile(work_dir + "/dft.log", dft_logger);
  }

  if (!do_dft_parse_ && (do_gwbse_ || do_esp_)) {
    // load the DFT data from serialized orbitals object
    std::string ORB_FILE =
        eqm_work_dir + "/molecules/" + frame_dir + "/" + orb_file;
    XTP_LOG(Log::error, pLog)
        << TimeStamp() << " Loading DFT data from " << ORB_FILE << std::flush;
    orbitals.ReadFromCpt(ORB_FILE);
  }

  if (do_gwbse_) {
    XTP_LOG(Log::error, pLog) << "Running GWBSE" << std::flush;
    try {
      GWBSE gwbse = GWBSE(orbitals);
      Logger gwbse_logger(votca::Log::current_level);
      gwbse_logger.setMultithreading(false);
      gwbse_logger.setPreface(Log::info, (format("\nGWBSE INF ...")).str());
      gwbse_logger.setPreface(Log::error, (format("\nGWBSE ERR ...")).str());
      gwbse_logger.setPreface(Log::warning, (format("\nGWBSE WAR ...")).str());
      gwbse_logger.setPreface(Log::debug, (format("\nGWBSE DBG ...")).str());
      gwbse.setLogger(&gwbse_logger);
      gwbse.Initialize(gwbse_options_);
      gwbse.Evaluate();
      gwbse.addoutput(segment_summary);
      WriteLoggerToFile(work_dir + "/gwbse.log", gwbse_logger);
    } catch (std::runtime_error& error) {
      std::string errormessage(error.what());
      SetJobToFailed(jres, pLog, "GWBSE:" + errormessage);
      return jres;
    }
  }

  if (do_esp_) {
    XTP_LOG(Log::error, pLog) << "Running ESPFIT" << std::flush;
    try {
      Esp2multipole esp2multipole = Esp2multipole(pLog);
      esp2multipole.Initialize(esp_options_);
      std::string ESPDIR =
          "MP_FILES/" + frame_dir + "/" + esp2multipole.GetStateString();
      StaticSegment seg2 = esp2multipole.Extractingcharges(orbitals);
      std::string mps_file = (format("%1%_%2%_%3%.mps") % segType % segId %
                              esp2multipole.GetStateString())
                                 .str();
      boost::filesystem::create_directories(ESPDIR);
      seg2.WriteMPS(ESPDIR + "/" + mps_file,
                    "Generated by eqm:" + esp2multipole.GetStateString());
      XTP_LOG(Log::error, pLog)
          << "Written charges to " << (ESPDIR + "/" + mps_file) << std::flush;
      segment_summary.add("partialcharges", (ESPDIR + "/" + mps_file));
    } catch (std::runtime_error& error) {
      std::string errormessage(error.what());
      SetJobToFailed(jres, pLog, "ESPFIT:" + errormessage);
      return jres;
    }
  }
  XTP_LOG(Log::error, pLog) << TimeStamp() << " Finished evaluating site "
                            << seg.getId() << std::flush;

  if (do_dft_parse_ || do_gwbse_) {
    XTP_LOG(Log::error, pLog) << "Saving data to " << orb_file << std::flush;
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
}  // namespace votca
