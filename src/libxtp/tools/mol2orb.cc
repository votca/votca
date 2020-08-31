#include "mol2orb.h"

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/xtp/molden.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

void Mol2Orb::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
  _moldenfile = _job_name + ".molden.input";
  _orbfile = _job_name + ".orb";

  _basisset_name = options.get(".basisset").as<std::string>();
  _aux_basisset_name = options.get(".auxbasisset").as<std::string>();
}

bool Mol2Orb::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Orbitals orbitals;

  Molden molden(_log);
  molden.setBasissetInfo(_basisset_name, _aux_basisset_name);
  molden.parseMoldenFile(_moldenfile, orbitals);

  // Save orbitals object
  XTP_LOG(Log::error, _log) << "Saving data to " << _orbfile << std::flush;
  orbitals.WriteToCpt(_orbfile);
  XTP_LOG(Log::error, _log) << "Done parsing\n" << std::flush;
  return true;
}

}  // namespace xtp
}  // namespace votca