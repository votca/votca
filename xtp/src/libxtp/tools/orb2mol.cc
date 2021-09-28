#include "orb2mol.h"

#include "votca/xtp/molden.h"
#include <boost/format.hpp>

namespace votca {
namespace xtp {

void Orb2Mol::ParseOptions(const tools::Property&) {

  moldenfile_ = job_name_ + ".molden.input";
  orbfile_ = job_name_ + ".orb";
  xyzfile_ = job_name_ + ".xyz";
}

bool Orb2Mol::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, log_) << "Loading data from " << orbfile_ << std::flush;
  orbitals.ReadFromCpt(orbfile_);

  XTP_LOG(Log::error, log_) << "Start parsing" << std::flush;

  Molden writer(log_);
  writer.WriteFile(moldenfile_, orbitals);

  XTP_LOG(Log::error, log_) << "Done parsing \n" << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca
