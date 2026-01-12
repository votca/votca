#ifndef VOTCA_XTP_EWALD_H
#define VOTCA_XTP_EWALD_H

#include "votca/xtp/IndexParser.h"
#include "votca/xtp/ewald/polarseg.h"
#include "votca/xtp/ewald/polartop.h"
#include <boost/format.hpp>
#include <boost/timer/timer.hpp>
#include <votca/xtp/ewald/ewald3d.h>
#include <votca/xtp/ewald/pewald3d.h>
#include <votca/xtp/ewald/xinductor.h>
#include <votca/xtp/ewald/xinteractor.h>
#include <votca/xtp/ewald/xjob.h>
#include <votca/xtp/ewald/xmapper.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/parallelxjobcalc.h>

#include <votca/xtp/ewald/polarbackground.h>

using boost::format;

namespace votca {
namespace xtp {

template <class EwaldMethod>
class Ewald : public ParallelXJobCalc<std::vector<Job>> {
 public:
  Ewald(){};
  ~Ewald(){};

  std::string Identify() const { return "ewald"; }
  void WriteJobFile(const Topology &top);
  void ReadJobFile(Topology &top);

  Job::JobResult EvalJob(const Topology &top, Job &job, QMThread &Thread);
  
  XJob ProcessInputString(Job &, const Topology &);

 protected:
  void ParseSpecificOptions(const tools::Property &user_options);

 private:
  tools::Property _options;
  // MULTIPOLES DEFINITION & MAPPING
  std::string _xml_file;
  std::string _mps_table;
  std::string _polar_bg_arch;
  std::string _mps_list_file;
  
  XMapper xmapper_;
  bool _pdb_check;
  bool _ptop_check;
  bool _use_mps_list;

  std::vector<QMState> states_;
  std::string which_segments_;
  Job createJob(const Segment &seg, const QMState &state, Index jobid) const;
};

template <class EwaldMethod>
void Ewald<EwaldMethod>::ParseSpecificOptions(const tools::Property &options) {
  _options = options;
  std::string key = ".multipoles";
  if (options.exists(key + ".polar_bg")) {
    _polar_bg_arch = options.get(key + ".polar_bg").as<std::string>();
  } else {
    _polar_bg_arch = "";
  }
  if (options.exists(key + ".pdb_check")) {
    _pdb_check = options.get(key + ".pdb_check").as<bool>();
  } else {
    _pdb_check = false;
  }
  if (options.exists(key + ".ptop_check")) {
    _ptop_check = options.get(key + ".ptop_check").as<bool>();
  } else {
    _ptop_check = false;
  }

  // for reading an optional list of segment specific MPS files
  if (options.exists(key + ".mps_table")) {
    _use_mps_list = options.get(key + ".mps_table").as<bool>();
  } else {
   _use_mps_list = false;
  }
  

  states_ = options.get("io_jobfile.states").as<std::vector<QMState>>();
  which_segments_ = options.get("io_jobfile.segments").as<std::string>();
}

template <class EwaldMethod>
Job Ewald<EwaldMethod>::createJob(const Segment &seg, const QMState &state,
                                  Index jobid) const {
  std::string marker = std::to_string(seg.getId()) + ":" + state.ToString();
  std::string tag = seg.getType() + "_" + marker;

  tools::Property Input;
  tools::Property &pInput = Input.add("input", "");
  pInput.add("site_energies", marker);
  pInput.add("state", state.ToString());
  pInput.add("segments", marker);
  Job job(jobid, tag, Input, Job::AVAILABLE);
  return job;
}

template <class EwaldMethod>
void Ewald<EwaldMethod>::ReadJobFile(Topology &top) {

  Index incomplete_jobs = 0;

  Eigen::Matrix<double, Eigen::Dynamic, 5> energies =
      Eigen::Matrix<double, Eigen::Dynamic, 5>::Zero(top.Segments().size(), 5);
  Eigen::Matrix<bool, Eigen::Dynamic, 5> found =
      Eigen::Matrix<bool, Eigen::Dynamic, 5>::Zero(top.Segments().size(), 5);

  tools::Property xml;
  xml.LoadFromXML(jobfile_);
  for (tools::Property *job : xml.Select("jobs.job")) {

    Index jobid = job->get("id").as<Index>();
    if (!job->exists("status")) {
      throw std::runtime_error(
          "Jobfile is malformed. <status> tag missing for job " +
          std::to_string(jobid));
    }
    if (job->get("status").as<std::string>() != "COMPLETE" ||
        !job->exists("output")) {
      incomplete_jobs++;
      continue;
    }

    std::vector<std::string> split =
        tools::Tokenizer(job->get("input.site_energies").as<std::string>(), ":")
            .ToVector();

    Index segid = std::stoi(split[0]);
    if (segid < 0 || segid >= Index(top.Segments().size())) {
      throw std::runtime_error("JobSegment id" + std::to_string(segid) +
                               " is not in topology for job " +
                               std::to_string(jobid));
    }
    QMState state;
    try {
      state.FromString(split[1]);
    } catch (std::runtime_error &e) {
      std::stringstream message;
      message << e.what() << " for job " << jobid;
      throw std::runtime_error(message.str());
    }
    double energy =
        job->get("output.summary.total").as<double>() * tools::conv::ev2hrt;
    if (found(segid, state.Type().Type()) != 0) {
      throw std::runtime_error("There are two entries in jobfile for segment " +
                               std::to_string(segid) +
                               " state:" + state.ToString());
    }

    energies(segid, state.Type().Type()) = energy;
    found(segid, state.Type().Type()) = true;
  }

  Eigen::Matrix<Index, 1, 5> found_states = found.colwise().count();
  std::cout << std::endl;
  for (Index i = 0; i < 5; i++) {
    if (found_states(i) > 0) {
      QMStateType type(static_cast<QMStateType::statetype>(i));
      std::cout << "Found " << found_states(i) << " states of type "
                << type.ToString() << std::endl;
    }
  }
  if (incomplete_jobs > 0) {
    std::cout << incomplete_jobs << " incomplete jobs found." << std::endl;
  }

  for (Segment &seg : top.Segments()) {
    Index segid = seg.getId();
    for (Index i = 0; i < 4; i++) {
      QMStateType type(static_cast<QMStateType::statetype>(i));
      if (found(segid, i) && found(segid, 4)) {
        double energy = energies(segid, i) - energies(segid, 4);
        seg.setEMpoles(type, energy);
      }
    }
  }

  return;
}

template <class EwaldMethod>
void Ewald<EwaldMethod>::WriteJobFile(const Topology &top) {

  std::cout << std::endl
            << "... ... Writing job file " << jobfile_ << std::flush;

  std::ofstream ofs;
  ofs.open(jobfile_, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("\nERROR: bad file handle: " + jobfile_);
  }

  std::vector<Index> segments_to_write;
  if (which_segments_ == "all") {
    for (Index i = 0; i < Index(top.Segments().size()); ++i) {
      segments_to_write.push_back(i);
    }
  } else {
    segments_to_write = IndexParser().CreateIndexVector(which_segments_);
  }

  ofs << "<jobs>" << std::endl;
  Index jobid = 0;
  for (Index segID : segments_to_write) {
    const Segment &seg = top.Segments()[segID];
    for (const QMState &state : states_) {
      Job job = createJob(seg, state, jobid);
      job.ToStream(ofs);
      jobid++;
    }
  }

  ofs << "</jobs>" << std::endl;
  ofs.close();
  std::cout << std::endl
            << "... ... In total " << jobid << " jobs" << std::flush;
}

template <class EwaldMethod>
XJob Ewald<EwaldMethod>::ProcessInputString(Job &job, const Topology &top) {

  std::vector<Segment *> qmSegs;
  std::vector<std::string> qmSegsState;

  // split <segments> tag at space
  std::vector<std::string> segments =
      tools::Tokenizer(job.getInput().get("segments").as<std::string>(), " ")
          .ToVector();

  for (auto segment : segments) {
    // now split each at :
    std::vector<std::string> split_segment =
        tools::Tokenizer(segment, ":").ToVector();
    Index segId = std::stoi(split_segment[0]);
    Segment seg = top.getSegment(segId);
    qmSegs.push_back(&seg);
    qmSegsState.push_back(segment);
  }

  return XJob(int(job.getId()), job.getTag(), qmSegs, qmSegsState, &top);
}

template <class EwaldMethod>
Job::JobResult Ewald<EwaldMethod>::EvalJob(const Topology &top, Job &job,
                                           QMThread &thread) {

  boost::timer::cpu_timer cpu_t;
  cpu_t.start();
  boost::timer::cpu_times t_in = cpu_t.elapsed();

  Logger &log = thread.getLogger();

  // CREATE XJOB FROM JOB INPUT std::string
  XJob xjob = this->ProcessInputString(job, top);
  XTP_LOG(Log::info, log) << "Created XJOB " << std::flush;

  // GENERATE POLAR TOPOLOGY (GENERATE VS LOAD IF PREPOLARIZED)
  if (_polar_bg_arch == "") {
    xmapper_.setLogger(&log);
    xmapper_.Gen_FGC_FGN_BGN(mapfile_, top, &xjob, _use_mps_list);
  } else {
    XTP_LOG(Log::info, log) << "Mps-Mapper: Generate FGC, load FGN BGN from '"
                            << _polar_bg_arch << "'" << std::flush;
    xmapper_.setLogger(&log);
    xmapper_.Gen_FGC_Load_FGN_BGN(mapfile_, top, &xjob, _polar_bg_arch, _use_mps_list);
  }

  // CALL THOLEWALD MAGIC
  EwaldMethod ewaldnd =
      EwaldMethod(&top, xjob.getPolarTop(), &_options, &thread.getLogger());

  if (_pdb_check) {
    ewaldnd.WriteDensitiesPDB(xjob.getTag() + ".densities.pdb");
  }

  ewaldnd.Evaluate();

  if (_ptop_check) {
    ewaldnd.WriteDensitiesPtop(xjob.getTag() + ".fg.ptop",
                               xjob.getTag() + ".bg.ptop");
  }

  // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
  tools::Property output = ewaldnd.GenerateOutputstring();
  Job::JobResult jres = Job::JobResult();
  jres.setOutput(output);
  jres.setStatus(Job::COMPLETE);

  if (!ewaldnd.Converged()) {
    jres.setStatus(Job::FAILED);
    jres.setError(ewaldnd.GenerateErrorstring());
    XTP_LOG(Log::error, log) << ewaldnd.GenerateErrorstring() << std::flush;
  }

  boost::timer::cpu_times t_out = cpu_t.elapsed();
  double t_run = (double(t_out.wall - t_in.wall)) / 1e9 / 60.;
  XTP_LOG(Log::info, log) << "Job runtime was " << t_run << " min"
                          << std::flush;

  return jres;
}

}  // namespace xtp
}  // namespace votca

#endif
