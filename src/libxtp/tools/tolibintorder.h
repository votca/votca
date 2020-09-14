
#pragma once
#ifndef VOTCA_XTP_TOLIBINTORDER_PRIVATE_H
#define VOTCA_XTP_TOLIBINTORDER_PRIVATE_H


// VOTCA includes
#include <votca/tools/constants.h>


// Libint include
#include <libint2.hpp>

// local VOTCA includes
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
  using real_t = libint2::scalar_type;
  using Matrix =
      Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Logger _log;
  Index _numofelectrons = 0;

  // Libint
  Matrix compute_1body_ints(
      const std::vector<libint2::Shell>& shells, libint2::Operator obtype,
      const std::vector<libint2::Atom>& atoms = std::vector<libint2::Atom>());
  size_t nbasis(const std::vector<libint2::Shell>& shells);
  size_t max_nprim(const std::vector<libint2::Shell>& shells);
  int max_l(const std::vector<libint2::Shell>& shells);
  std::vector<size_t> map_shell_to_basis_function(
      const std::vector<libint2::Shell>& shells);

  std::vector<libint2::Shell> make_libint_basis(const AOBasis& aobasis);

  // VOTCA
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOLIBINTORDER_PRIVATE_H