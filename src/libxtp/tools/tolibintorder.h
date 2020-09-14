
#pragma once
#ifndef VOTCA_XTP_TOLIBINTORDER_PRIVATE_H
#define VOTCA_XTP_TOLIBINTORDER_PRIVATE_H

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
  ToLibintOrder() = default;

  ~ToLibintOrder() final = default;

  std::string Identify() final { return "tolibintorder"; }

  void Initialize(const tools::Property& user_options) final;
  bool Evaluate() final;

 private:
  // clang-format off
 std::array<Index,25> _libint_multipliers={ 
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            1,1,1,1,1,1,1, //f 
            1,1,1,1,1,1,1,1,1 //g
            };
  std::array<Index, 25> _libint_reorder={ 
            0, //s
            -1,0,1, //p
            -2,-1,0,1,2, //d
            -3,-2,-1,0,1,2,3, //f 
            -4,-3,-2,-1,0,1,2,3,4 //g
            };
  // clang-format on
  Logger _log;
  Index _numofelectrons = 0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOLIBINTORDER_PRIVATE_H