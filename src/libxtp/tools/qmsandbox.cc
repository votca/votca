#include "qmsandbox.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/threecenter.h"
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

  BasisSet auxbs;
  auxbs.Load("aux-def2-svp");

  AOBasis auxbasis;
  auxbasis.Fill(auxbs, atoms);

  std::vector<libint2::Shell> shells = basis.GenerateLibintBasis();
  std::vector<libint2::Shell> auxShells = auxbasis.GenerateLibintBasis();

  std::copy(std::begin(shells), std::end(shells),
            std::ostream_iterator<libint2::Shell>(std::cout, "\n"));

  std::cout << "*************************************************************\n"
            << std::endl;

  std::copy(std::begin(auxShells), std::end(auxShells),
            std::ostream_iterator<libint2::Shell>(std::cout, "\n"));

  std::cout << "*************************************************************\n"
            << std::endl;

  TCMatrix_dft three;
  three.Fill(auxbasis, basis);

  TCMatrix_dft three2;
  three2.Fill2(auxbasis, basis);

  std::cout << "SIZE " << three.size() << std::endl;

  for (int i = 0; i < three.size(); i++) {
    std::cout << three[i] << std::endl;
    std::cout << three2[i] << std::endl;
    std::cout << "*************************************************************"
              << std::endl;
  }

  return true;
}  // namespace xtp

}  // namespace xtp
}  // namespace votca