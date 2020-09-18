#include "qmsandbox.h"

namespace votca {
namespace xtp {

void QMSandbox::Initialize(const tools::Property& user_options) {

  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
}

bool QMSandbox::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  QMMolecule atoms("", 0);
  atoms.LoadFromFile(_job_name + ".xyz");

  BasisSet bs;
  bs.Load("def2-tzvp");
  AOBasis basis;
  basis.Fill(bs, atoms);

  // Compute Overlap with VOTCA
  AOOverlap dftAOoverlap;
  dftAOoverlap.Fill(basis);

  XTP_LOG(Log::error, _log) << "\n\tOverlap Integrals VOTCA:\n";
  XTP_LOG(Log::error, _log) << dftAOoverlap.Matrix() << std::endl;

  return true;
}

}  // namespace xtp
}  // namespace votca