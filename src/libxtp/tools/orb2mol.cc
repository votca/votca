#include "orb2mol.h"

#include "votca/xtp/moldenwriter.h"
#include <boost/format.hpp>

namespace votca {
namespace xtp {

void Orb2Mol::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
  _moldenfile = _job_name + ".molden.input";
  _orbfile = _job_name + ".orb";
  _xyzfile = _job_name + ".xyz";
}

bool Orb2Mol::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, _log) << "Loading data from " << _orbfile << std::flush;
  orbitals.ReadFromCpt(_orbfile);

  XTP_LOG(Log::error, _log) << "Start parsing" << std::flush;

  MoldenWriter writer(_log);
  writer.WriteFile(_moldenfile, orbitals);

  XTP_LOG(Log::error, _log) << "Done parsing \n" << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca
