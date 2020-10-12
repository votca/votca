#include "qmsandbox.h"
#include "votca/xtp/orbreorder.h"
#include <votca/tools/eigenio_matrixmarket.h>

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

  OPENMP::setMaxThreads(_nThreads);

  QMMolecule atoms("", 0);
  atoms.LoadFromFile(_job_name + ".xyz");

  BasisSet bs;
  bs.Load("def2-svp");
  AOBasis basis;
  basis.Fill(bs, atoms);

  std::vector<libint2::Shell> shells = basis.GenerateLibintBasis();

  for (auto& shell : basis) {
    std::cout << shell << std::endl;
  }

  std::copy(std::begin(shells), std::end(shells),
            std::ostream_iterator<libint2::Shell>(std::cout, "\n"));

  AOOverlap aooverlap;
  AOCoulomb aocoulomb;
  AOKinetic aokinetic;
  AODipole aodipole;

  for (int i = 0; i < 100; i++) {
    aooverlap.Fill(basis);
  }
  for (int i = 0; i < 100; i++) {
    aocoulomb.Fill(basis);
  }
  for (int i = 0; i < 100; i++) {
    aokinetic.Fill(basis);
  }
  for (int i = 0; i < 100; i++) {
    aodipole.Fill(basis);
  }

  return true;
}  // namespace xtp

}  // namespace xtp
}  // namespace votca