
#pragma once
#ifndef VOTCA_XTP_TOLIBINTORDER_H
#define VOTCA_XTP_TOLIBINTORDER_H

// VOTCA includes
#include <votca/tools/constants.h>

// Libint include
#include <libint2.hpp>

// local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

class ToLibintOrder final : public QMTool {
 public:
  ToLibintOrder() { _libint_multipliers.fill(1); }

  ~ToLibintOrder() final = default;

  std::string Identify() final { return "tolibintorder"; }

  void Initialize(const tools::Property& user_options) final;
  bool Evaluate() final;

 private:
  // clang-format off
 std::array<Index,49> _libint_multipliers;
  std::array<Index, 49> _libint_reorder={
      0,                             // s
      0, -1, 1,                      // p
      0, -1, 1, -2, 2,               // d
      0, -1, 1, -2, 2, -3, 3,        // f
      0, -1, 1, -2, 2, -3, 3, -4, 4,  // g
      0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5,  // h
      0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6  // i
  };

  // clang-format on
  Logger _log;
  std::string _xyz_file;
  std::string _basis_set;
  std::string _matrix_file;
  bool _aux_basis;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOLIBINTORDER_H