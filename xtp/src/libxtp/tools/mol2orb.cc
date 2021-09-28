#include "mol2orb.h"

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/xtp/molden.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

void Mol2Orb::ParseOptions(const tools::Property& options) {

  moldenfile_ = job_name_ + ".molden.input";
  orbfile_ = job_name_ + ".orb";

  basisset_name_ = options.get(".basisset").as<std::string>();
  aux_basisset_name_ = options.get(".auxbasisset").as<std::string>();
}

bool Mol2Orb::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  Orbitals orbitals;

  Molden molden(log_);
  molden.setBasissetInfo(basisset_name_, aux_basisset_name_);
  molden.parseMoldenFile(moldenfile_, orbitals);

  // Save orbitals object
  XTP_LOG(Log::error, log_) << "Saving data to " << orbfile_ << std::flush;
  orbitals.WriteToCpt(orbfile_);
  XTP_LOG(Log::error, log_) << "Done parsing\n" << std::flush;
  return true;
}

}  // namespace xtp
}  // namespace votca