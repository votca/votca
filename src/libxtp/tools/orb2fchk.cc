#include "orb2fchk.h"

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/xtp/gaussianwriter.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

void Orb2Fchk::ParseOptions(const tools::Property& options) {

  _fchk_filename = _job_name + ".fchk";
  _orbfile = _job_name + ".orb";
}

bool Orb2Fchk::Run() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, _log) << "Loading data from " << _orbfile << std::flush;
  orbitals.ReadFromCpt(_orbfile);

  XTP_LOG(Log::error, _log) << "Start parsing" << std::flush;

  GaussianWriter writer(_log);
  writer.WriteFile(_fchk_filename, orbitals);

  XTP_LOG(Log::error, _log) << "Done parsing \n" << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca