#ifndef VOTCA_XTP_EWDBGPOL_H
#define VOTCA_XTP_EWDBGPOL_H

#include <boost/filesystem.hpp>
#include <votca/tools/tokenizer.h>
#include <votca/xtp/ewald/polarbackground.h>
#include <votca/xtp/ewald/qmthread.h>
// #include <votca/xtp/ewald/xmapper.h>
#include "votca/xtp/backgroundregion.h"
#include "votca/xtp/ewald/polarseg.h"
#include "votca/xtp/ewald/polartop.h"
#include "votca/xtp/segmentmapper.h"
#include <votca/xtp/job.h>
#include <votca/xtp/jobtopology.h>
#include <votca/xtp/qmcalculator.h>
#include <votca/tools/getline.h>

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
  std::string _xml_file;
  bool _use_mps_table;
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

  // std::string key = "options.ewdbgpol";
  std::string key = "";
  if (opt.exists(key + ".multipoles")) {
    _xml_file = opt.get(key + ".multipoles").as<std::string>();
  } else {
    std::cout << std::endl;
    throw std::runtime_error("No multipole mapping file provided");
  }

  key = ".control";
  // CONTROL
  if (opt.exists(key + ".mps_table")) {
    _use_mps_table = opt.get(key + ".mps_table").as<bool>();
  } else {
    _use_mps_table = false;
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

  // use new MAPPER to get a list of mapped PolarSegments (neutral only for now)
  BackgroundRegion BGN(0, log);
  PolarMapper polmap(log);
  polmap.LoadMappingFile(_xml_file);
  Index seg_index = 0;

  // if mps_table file is given:
  // - read it
  // use the filenames in polmap.map
  if (_use_mps_table) {
    // MPS mapping is performed with a unique MPS file for each segment.
    // filenames are expected to be "MP_FILES/segmentname_segmentid_state.mps"
    // e.g., MP_FILES/edot_1_n.mps
    for (auto segment : top.Segments()) {    
      std::string mpsfile_n = "MP_FILES/" + segment.getType() + "_" + std::to_string(segment.getId()) + "_n.mps";
      std::cout << mpsfile_n << std::endl;
      PolarSegment mol = polmap.map(segment, mpsfile_n);
      BGN.push_back(mol);
      seg_index++;
    }

      exit(0);

  } else {
    for (auto segment : top.Segments()) {
      PolarSegment mol = polmap.map(segment, SegId(seg_index, "n"));
      BGN.push_back(mol);
      seg_index++;
    }
  }

  // Convert this to old PolarTop
  PolarTop ptop(&top);

  // DECLARE TARGET CONTAINERS
  std::vector<PolarSeg *> bgN;
  std::vector<Segment *> segs_bgN;

  bgN.reserve(BGN.size());
  // PARTITION SEGMENTS ONTO BACKGROUND + FOREGROUND
  segs_bgN.reserve(BGN.size());
  for (auto segment : top.Segments()) {
    segs_bgN.push_back(&segment);
  }

  // segments in BGN are NEW POLARSEGMENTS
  int state = 0;  // only neutral background segments
  for (auto segment : BGN) {
    // get all NEW PolarSites of this segment and convert them to OLD PolarSites
    std::vector<APolarSite *> psites;
    psites.reserve(segment.size());
    for (auto site : segment) {
      APolarSite *psite = new APolarSite();
      psite->ConvertFromPolarSite(site, state);
      psite->Charge(state);  // set state of this sites 0: ground state
      psites.push_back(psite);
    }

    // now make an OLD PolarSeg from the new PolarSegment
    PolarSeg *new_pseg = new PolarSeg(int(segment.getId()), psites);
    // std::shared_ptr<PolarSeg> new_pseg( new  PolarSeg(int(segment.getId()),
    // psites));

    bgN.push_back(new_pseg);
  }

  // PROPAGATE SHELLS TO POLAR TOPOLOGY
  ptop.setBGN(bgN);
  ptop.setSegsBGN(segs_bgN);

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
  EWD::PolarBackground pbg(&top, &ptop, _options, &log);
  pbg.Polarize(int(nThreads_));

  // SAVE POLARIZATION STATE
  if (pbg.HasConverged()) {
    XTP_LOG(Log::info, log) << "Save polarization state" << std::flush;
    ptop.SaveToDrive("bgp_main.ptop");
    ptop.PrintPDB("bgp_main.pdb");
  }

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
