#ifndef VOTCA_XTP_EWALD_H
#define VOTCA_XTP_EWALD_H

#include <votca/xtp/parallelxjobcalc.h>
#include <votca/xtp/ewald/xinductor.h>
#include <votca/xtp/ewald/xinteractor.h>
#include <votca/xtp/ewald/xjob.h>
#include <boost/format.hpp>
#include <boost/timer/timer.hpp>
#include <votca/xtp/ewald/ewald3d.h>
#include <votca/xtp/ewald/pewald3d.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/ewald/xmapper.h>
#include "votca/xtp/ewald/polarseg.h"
#include "votca/xtp/ewald/polartop.h"

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
  void Initialize(tools::Property *);
  void WriteJobFile(const Topology &top);
  void ReadJobFile(Topology &top);

  void PreProcess(Topology *top);
  Job::JobResult EvalJob(const Topology &top, Job &job, QMThread &Thread);
  void PostProcess(Topology *top) { ; }

  XJob ProcessInputString(Job &, const Topology &, QMThread &);

 protected:
  void ParseSpecificOptions(const tools::Property &user_options);

 private:
  tools::Property _options;
  // MULTIPOLES DEFINITION & MAPPING
  std::string _xml_file;
  std::string _mps_table;
  std::string _polar_bg_arch;
  XMapper                        xmapper_;
  bool _pdb_check;
  bool _ptop_check;
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

  /* print_regions_pdb_ = options.get(".print_regions_pdb").as<bool>();
  max_iterations_ = options.get(".max_iterations").as<Index>();
  regions_def_.second = options.get(".regions");
  regions_def_.first = mapfile_;
  use_gs_for_ex_ = options.get("io_jobfile.use_gs_for_ex").as<bool>();

  states_ = options.get("io_jobfile.states").as<std::vector<QMState>>();
  which_segments_ = options.get("io_jobfile.segments").as<std::string>();

  bool groundstate_found = std::any_of(
      states_.begin(), states_.end(),
      [](const QMState& state) { return state.Type() == QMStateType::Gstate; });
  if (!groundstate_found) {
    states_.push_back(QMState("n"));
  }*/
}

template <class EwaldMethod>
void Ewald<EwaldMethod>::Initialize(tools::Property *opt) {

  // NOTE These options are passed on to <EwaldMethod> in ::EvalJob
  //_options = opt;
  //_maverick = (_nThreads == 1) ? true : false;

  /* std::cout << std::endl
       << "... ... Initialized with " << nThreads_ << " threads. "
       << std::endl;

  std::string key = ".jobcontrol";
      if ( opt->exists(key+".job_file")) {
          jobfile_ = opt->get(key+".job_file").as<std::string>();
      }
      else {
          std::cout << std::endl;
          throw std::runtime_error("Job-file not set. Abort.");
      }

  key = "options.ewald.multipoles";
      if (opt->exists(key+".mapping")) {
          _xml_file = opt->get(key+".mapping").as< std::string >();
      }
      else {
          std::cout << std::endl;
          throw std::runtime_error("Multipole mapping file not set. Abort.");
      }
      if ( opt->exists(key+".mps_table")) {
          _mps_table = opt->get(key+".mps_table").as<std::string>();
      }
      else {
          std::cout << std::endl;
          throw std::runtime_error("Background mps table not set. Abort.");
      }

      */

  return;
}

template <class EwaldMethod>
void Ewald<EwaldMethod>::PreProcess(Topology *top) {
  // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
  std::cout << std::endl << "... ... Initialize MPS-mapper: " << std::flush;
  //_mps_mapper.GenerateMap(_xml_file, _mps_table, top);
  return;
}

template <class EwaldMethod>
void Ewald<EwaldMethod>::ReadJobFile(Topology &top) {
  assert(false && "<::ReadJobFile> NOT IMPLEMENTED");
  return;
}

template <class EwaldMethod>
void Ewald<EwaldMethod>::WriteJobFile(const Topology &top) {

  // SET UP FILE STREAM
  std::ofstream ofs;
  std::string jobFile = jobfile_;
  ofs.open(jobFile.c_str(), std::ofstream::out);
  if (!ofs.is_open()) throw std::runtime_error("Bad file handle: " + jobFile);

  ofs << "<jobs>" << std::endl;

  int jobCount = -1;
  // std::vector<Segment*>::iterator sit1;

  // DEFINE PAIR CHARGE STATES
  std::vector<std::string> states;
  std::vector<std::string>::iterator vit;
  states.push_back("n");
  states.push_back("e");
  states.push_back("h");

  // CREATE JOBS FOR ALL SEGMENTS AND STATES
  for (auto &seg1 : top.Segments()) {

    int id1 = seg1.getId();
    std::string name1 = seg1.getType();

    for (vit = states.begin(); vit != states.end(); ++vit) {
      int id = ++jobCount;
      std::string s1 = *vit;
      std::string tag = (format("%1$d:%3$s:%2$s") % id1 % s1 % name1).str();

      tools::Property input;
      tools::Property &out = input.add("input", "");
      tools::Property *next = NULL;
      next = &out.add("segment", "");
      next->add("id", (format("%1$d") % id1).str());
      next->add("type", (format("%1$s") % name1).str());
      next->add("mps", (format("MP_FILES/%1$s_%2$s.mps") % name1 % s1).str());

      Job job(id, tag, input, Job::AVAILABLE);
      job.ToStream(ofs);
    }
  }

  // CLOSE STREAM
  ofs << "</jobs>" << std::endl;
  ofs.close();
}

template <class EwaldMethod>
XJob Ewald<EwaldMethod>::ProcessInputString(Job &job, const Topology &top,
                                            QMThread &thread) {

  // Input std::string looks like this:
  // <id1>:<name1>:<mpsfile1> <id2>:<name2>: ... ... ...
  XTP_LOG(Log::error, thread.getLogger())
      << "EWALD::ProcessInput String currently only for single segment! "
      << std::flush;
  std::vector<Segment *> qmSegs;
  std::vector<std::string> qmSegMps;

  Index segId = job.getInput().get("segment.id").as<Index>();
  std::string segName = job.getInput().get("segment.type").as<std::string>();

  Segment seg = top.getSegment(segId);
  if (seg.getType() != segName) {
    XTP_LOG(Log::error, thread.getLogger())
        << "ERROR: Seg " << segId << ":" << seg.getType() << " "
        << " maltagged as " << segName << ". Skip job ..." << std::flush;
    throw std::runtime_error("Input does not match topology.");
  }

  std::string mpsFile = job.getInput().get("segment.mps").as<std::string>();
  qmSegs.push_back(&seg);
  qmSegMps.push_back(mpsFile);

  return XJob(job.getId(), job.getTag(), qmSegs, qmSegMps, &top);
}

template <class EwaldMethod>
Job::JobResult Ewald<EwaldMethod>::EvalJob(const Topology &top, Job &job,
                                           QMThread &thread) {

  boost::timer::cpu_timer cpu_t;
  cpu_t.start();
  boost::timer::cpu_times t_in = cpu_t.elapsed();

  Logger &log = thread.getLogger();
  XTP_LOG(Log::info, log)
      //<< "Job input = " << job.getInput().as<std::string>() << std::flush;
      << "Job input = " << job.getInput() << std::flush;

  // CREATE XJOB FROM JOB INPUT std::string
  XJob xjob = this->ProcessInputString(job, top, thread);
  XTP_LOG(Log::info, log) << "Created XJOB " << std::flush;

  // GENERATE POLAR TOPOLOGY (GENERATE VS LOAD IF PREPOLARIZED)
  // use new MAPPER to get a list of mapped PolarSegments (neutral only for now)

  // Convert this to old PolarTop
  //Topology new_top = top;
  //PolarTop ptop(&new_top);
  if (_polar_bg_arch == "") {
        xmapper_.setLogger(&log);
        xmapper_.Gen_FGC_FGN_BGN(mapfile_, top, &xjob);
  } else {
    XTP_LOG(Log::info, log) << "LOADING FROM ARCHIVE TODO!!!!! '"
                            << _polar_bg_arch << "'" << std::flush;
    XTP_LOG(Log::info, log) << "Mps-Mapper: Generate FGC, load FGN BGN from '"
                            << _polar_bg_arch << "'" << std::flush;
    //_mps_mapper.Gen_FGC_Load_FGN_BGN(top, &xjob, _polar_bg_arch, thread);
  }

  // CALL THOLEWALD MAGIC
  XTP_LOG(Log::info, log) << "Trying to construct ewald object" << std::flush;

  //std::cout << ptop.FGN().size() << "\n" << std::endl;
  //std::cout << ptop.FGC().size() << "\n" << std::endl;
  //std::cout << ptop.BGN().size() << "\n" << std::endl;

  EwaldMethod ewaldnd =
      EwaldMethod(&top, xjob.getPolarTop(), &_options, &thread.getLogger());
  XTP_LOG(Log::info, log) << "success" << std::flush;

  if (_pdb_check) ewaldnd.WriteDensitiesPDB(xjob.getTag() + ".densities.pdb");
  ewaldnd.Evaluate();
  if (_ptop_check)
    ewaldnd.WriteDensitiesPtop(xjob.getTag() + ".fg.ptop",
                               xjob.getTag() + ".mg.ptop",
                               xjob.getTag() + ".bg.ptop");

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
  double t_run = (t_out.wall - t_in.wall) / 1e9 / 60.;
  XTP_LOG(Log::info, log) << "Job runtime was " << t_run << " min"
                          << std::flush;

  return jres;
}

}  // namespace xtp
}  // namespace votca

#endif
