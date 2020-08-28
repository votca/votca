
#pragma once
#ifndef VOTCA_XTP_MOL2ORB_PRIVATE_H
#define VOTCA_XTP_MOL2ORB_PRIVATE_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

class Mol2Orb final : public QMTool {
 public:
  Mol2Orb() = default;

  ~Mol2Orb() final = default;

  std::string Identify() final { return "mol2orb"; }

  void Initialize(const tools::Property& user_options) final;
  bool Evaluate() final;

 private:
  std::string _moldenfile;
  std::string _orbfile;
  std::string _xyzfile;
  std::string _basisset_name;
  std::string _aux_basisset_name;
  AOBasis _basis;
  AOBasis _auxbasis;
  Logger _log;

  std::string readAtoms(QMMolecule& mol, const std::string& units,
                        std::ifstream& input_file) const;
  std::string readMOs(Orbitals& orbitals, std::ifstream& input_file);
  void addBasissetInfo(Orbitals& orbitals);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOL2ORB_PRIVATE_H