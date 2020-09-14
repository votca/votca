#include "tolibintorder.h"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/orbreorder.h"

namespace votca {
namespace xtp {

void ToLibintOrder::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);

  _xyz_file =
      options.ifExistsReturnElseReturnDefault<std::string>("xyz", _job_name);

  _basis_set = options.ifExistsReturnElseReturnDefault<std::string>("basisset",
                                                                    _job_name);

  _matrix_file =
      options.ifExistsReturnElseReturnDefault<std::string>("matrix", _job_name);
}

bool ToLibintOrder::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Eigen::MatrixXd inputMatrix =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(_matrix_file);

  QMMolecule atoms("", 0);
  atoms.LoadFromFile(_xyz_file);

  BasisSet bs;
  bs.Load(_basis_set);
  AOBasis basis;
  basis.Fill(bs, atoms);

  OrbReorder reorder(_libint_reorder, _libint_multipliers);

  reorder.reorderOperator(inputMatrix, basis);

  votca::tools::EigenIO_MatrixMarket::WriteMatrix("output.mm", inputMatrix);

  XTP_LOG(Log::error, _log) << _job_name + ".mm" << std::endl;

  return true;
}

}  // namespace xtp
}  // namespace votca