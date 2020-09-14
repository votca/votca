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
}

bool ToLibintOrder::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Eigen::MatrixXd inputMatrix =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(_job_name);

  OrbReorder reorder(_libint_reorder, _libint_multipliers);

  reorder.reorderRowsAndCols(inputMatrix, inputBasis);

  votca::tools::EigenIO_MatrixMarket::WriteMatrix("output.mm", inputMatrix);

  XTP_LOG(Log::error, _log) << _job_name + ".mm" << std::endl;

  return true;
}

}  // namespace xtp
}  // namespace votca