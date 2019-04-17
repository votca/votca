/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_XTP_GWBSE_H
#define VOTCA_XTP_GWBSE_H
#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/gw.h>
#include <votca/xtp/logger.h>

#include "bse.h"

namespace votca {
namespace xtp {
class Orbitals;
class AOBasis;
/**
 * \brief Electronic excitations from GW-BSE
 *
 * Evaluates electronic excitations in molecular systems based on
 * many-body Green's functions theory within the GW approximation and
 * the Bethe-Salpeter equation. Requires molecular orbitals
 *
 *  B. Baumeier, Y. Ma, D. Andrienko, M. Rohlfing
 *  J. Chem. Theory Comput. 8, 997-1002 (2012)
 *
 *  B. Baumeier, D. Andrienko, M. Rohlfing
 *  J. Chem. Theory Comput. 8, 2790-2795 (2012)
 *
 */

class GWBSE {
 public:
  GWBSE(Orbitals& orbitals) : _orbitals(orbitals){};

  void Initialize(tools::Property& options);

  std::string Identify() { return "gwbse"; }

  void setLogger(Logger* pLog) { _pLog = pLog; }

  bool Evaluate();

  void addoutput(tools::Property& summary);

 private:
  Eigen::MatrixXd CalculateVXC(const AOBasis& dftbasis);
  int CountCoreLevels();
  Logger* _pLog;
  Orbitals& _orbitals;

  // program tasks
  bool _do_qp_diag;
  bool _do_bse_diag;
  bool _do_bse_singlets;
  bool _do_bse_triplets;

  // storage tasks
  bool _store_qp_pert;
  bool _store_qp_diag;
  bool _store_bse_singlets;
  bool _store_bse_triplets;
  bool _store_eh_interaction;

  // options for own Vxc calculation
  bool _doVxc;
  std::string _functional;
  std::string _grid;

  int _openmp_threads;

  // fragment definitions
  int _fragA;

  // BSE variant

  GW::options _gwopt;
  BSE::options _bseopt;

  // basis sets
  std::string _auxbasis_name;
  std::string _dftbasis_name;

  std::vector<QMFragment<BSE_Population> > _triplets;
  std::vector<QMFragment<BSE_Population> > _singlets;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GWBSE_H
