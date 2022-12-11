#ifndef VOTCA_XTP_EWDBGPOL_H
#define VOTCA_XTP_EWDBGPOL_H

#include <boost/filesystem.hpp>
#include <votca/tools/tokenizer.h>
#include <votca/xtp/ewald/polarbackground.h>
#include <votca/xtp/ewald/qmthread.h>
// #include <votca/xtp/ewald/xmapper.h>
#include "votca/xtp/backgroundregion.h"
#include "votca/xtp/segmentmapper.h"
#include <votca/xtp/job.h>
#include <votca/xtp/jobtopology.h>
#include <votca/xtp/qmcalculator.h>

namespace votca {
namespace xtp {

class EwaldBgPolarizer final : public QMCalculator {
 public:
  std::string Identify() const { return "ewdbgpol"; }
  bool WriteToStateFile() const { return false; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Evaluate(Topology &top);

 private:
  tools::Property _options;
  std::string _mps_table;
  std::string _xml_file;
  // XMpsMap _mps_mapper;
  bool _pdb_check;

  std::string _ptop_file;
  bool _do_restart;
};

void EwaldBgPolarizer::ParseOptions(const tools::Property &opt) {

  _options = opt;
  //_maverick = (_nThreads == 1) ? true : false;
  std::cout << std::endl
            << "... ... Initialized with " << nThreads_ << " threads. "
            << std::flush;

  std::string key = "options.ewdbgpol";
  if (opt.exists(key + ".multipoles")) {
    _xml_file = opt.get(key + ".multipoles").as<std::string>();
  } else {
    std::cout << std::endl;
    throw std::runtime_error("No multipole mapping file provided");
  }

  key = "options.ewdbgpol.control";
  // CONTROL
  if (opt.exists(key + ".mps_table")) {
    _mps_table = opt.get(key + ".mps_table").as<std::string>();
  } else {
    _mps_table = opt.get(key + ".emp_file").as<std::string>();
  }
  if (opt.exists(key + ".pdb_check")) {
    _pdb_check = opt.get(key + ".pdb_check").as<bool>();
  } else {
    _pdb_check = false;
  }
  // RESTART OPTIONS
  if (opt.exists(key + ".restart_from")) {
    _ptop_file = opt.get(key + ".restart_from").as<std::string>();
    if (boost::filesystem::exists(_ptop_file))
      _do_restart = true;
    else {
      _do_restart = false;
      std::cout << std::endl
                << "... ... File " << _ptop_file
                << " not found. Ignoring Restart_from option" << std::flush;
    }
  } else {
    _ptop_file = "";
    _do_restart = false;
  }

  return;
}

bool EwaldBgPolarizer::Evaluate(Topology &top) {

  QMThread master;
  master.getLogger()->setReportLevel(Log::debug);
  master.getLogger()->setMultithreading(true);
  master.getLogger()->setPreface(Log::info, "\nMST INF");
  master.getLogger()->setPreface(Log::error, "\nMST ERR");
  master.getLogger()->setPreface(Log::warning, "\nMST WAR");
  master.getLogger()->setPreface(Log::debug, "\nMST DBG");
  Logger &log = *master.getLogger();

  // GENERATE BACKGROUND (= periodic bg, with empty foreground)
  std::cout << std::endl << "... ... Initialize MPS-mapper: " << std::flush;

  // * THIS NEEDS TO BE ADAPTED TO NEW JOBTOPOLOGY instead of polartopology
  // first define a job
  // tools::Property jobproperty;
  // jobproperty.add("id",0);
  // jobproperty.add("tag", "background");
  // jobproperty.add("input","background");
  // jobproperty.add("status","AVAILABLE");
  // Job job(jobproperty);
  // JobTopology jobtop = JobTopology(job, log, "EWDBACKGROUND");
  BackgroundRegion BGN(0, log);

  PolarMapper polmap(log);
  polmap.LoadMappingFile(_xml_file);
  Index seg_index = 0;
  for (auto segment : top.Segments()) {

    PolarSegment mol = polmap.map(segment, SegId(seg_index, "n"));
    BGN.push_back(mol);
    seg_index++;
  }
  // for (const SegId& seg_index : seg_ids) {
  //  const Segment& segment = top.getSegment(seg_index.Id());

  // PolarSegment mol = polmap.map(segment, seg_index);

  // ShiftPBC(top, center, mol);
  // mol.setType("mm" + std::to_string(id));
  // polarregion->push_back(mol);

  /*PolarTop ptop(&top);
  if (_do_restart) {
    ptop.LoadFromDrive(_ptop_file);
  } else {
    _mps_mapper.GenerateMap(_xml_file, _mps_table, &top);
    _mps_mapper.Gen_BGN(&top, &ptop, &master);
  }
*/

  if (_pdb_check) {
    csg::PDBWriter mpwriter;
    mpwriter.Open("ewdbgpol.pdb", false);
    mpwriter.WriteHeader("Background");
    mpwriter.WriteBox(top.getBox() * tools::conv::bohr2ang);
    BGN.WritePDB(mpwriter);
  }

  // POLARIZE SYSTEM
  EWD::PolarBackground pbg(&top, &BGN, &_options, &log);
  // pbg.Polarize(nThreads_);

  // SAVE POLARIZATION STATE
  // if (pbg.HasConverged()) {
  //  XTP_LOG(Log::info, log) << "Save polarization state" << std::flush;
  //  ptop.SaveToDrive("bgp_main.ptop");
  //  ptop.PrintPDB("bgp_main.pdb");
  //}

  //    // LOAD POLARIZATION STATE
  //    CTP_LOG(logINFO,log) << "Load polarization state" << flush;
  //    PolarTop ptop2;
  //    ptop2.LoadFromDrive("bgp_main.ptop");
  //    ptop2.PrintPDB("bgp_check.pdb");
  //    ptop2.SaveToDrive("bgp_check.ptop");

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif
