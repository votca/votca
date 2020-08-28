
#pragma once
#ifndef VOTCA_XTP_ORB2MOL_PRIVATE_H
#define VOTCA_XTP_ORB2MOL_PRIVATE_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <fstream>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

class Orb2Mol final : public QMTool {
 public:
  Orb2Mol() = default;

  ~Orb2Mol() final = default;

  std::string Identify() final { return "orb2mol"; }

  void Initialize(const tools::Property& user_options) final;
  bool Evaluate() final;

 private:
  // clang-format off
  std::array<Index,25> _multipliers={
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            1,1,1,1,1,-1,-1, //f 
            1,1,1,1,1,-1,-1,-1,-1 //g
            };

  OrbTranspositions _transpositions { 
    std::vector<std::array<Index, 2>> {}, //s
    std::vector<std::array<Index, 2>> {   //p
      {0, 2}
    }, 
    std::vector<std::array<Index, 2>> {   //d
      {1, 2},
      {3, 4}
      }, 
    std::vector<std::array<Index, 2>> {   //f
      {1, 2},  
      {3, 4},
      {5, 6}
    }, 
    std::vector<std::array<Index, 2>> {   //g
      {1, 2},
      {3, 4},
      {5, 6},
      {7, 8}
    }
  };
  std::string _moldenfile;
  std::string _orbfile;
  std::string _xyzfile;
  std::string _basisset_name;
  std::string _aux_basisset_name;
  AOBasis _basis;
  BasisSet _bs;
  Logger _log;

 void writeAtoms(Orbitals& orbitals, std::ofstream& outFile);
 void writeMOs(Orbitals& orbitals, std::ofstream& outFile);
 void writeBasisSet(Orbitals& orbitals, std::ofstream& outFile);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORB2MOL_PRIVATE_H 
