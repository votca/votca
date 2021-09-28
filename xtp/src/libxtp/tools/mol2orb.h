
#pragma once
#ifndef VOTCA_XTP_MOL2ORB_H
#define VOTCA_XTP_MOL2ORB_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

class Mol2Orb final : public QMTool {
 public:
  Mol2Orb() = default;

  ~Mol2Orb() = default;

  std::string Identify() { return "mol2orb"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  std::string _moldenfile;
  std::string _orbfile;
  std::string _basisset_name;
  std::string _aux_basisset_name;
  Logger _log;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOL2ORB_H