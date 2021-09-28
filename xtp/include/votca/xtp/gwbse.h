/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#pragma once
#ifndef VOTCA_XTP_GWBSE_H
#define VOTCA_XTP_GWBSE_H

// Standard includes
#include <fstream>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "bse.h"
#include "eigen.h"
#include "gw.h"
#include "logger.h"
#include "qmfragment.h"

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
  GWBSE(Orbitals& orbitals) : orbitals_(orbitals){};

  void Initialize(tools::Property& options);

  std::string Identify() { return "gwbse"; }

  void setLogger(Logger* pLog) { pLog_ = pLog; }

  bool Evaluate();

  void addoutput(tools::Property& summary);

 private:
  Eigen::MatrixXd CalculateVXC(const AOBasis& dftbasis);
  Index CountCoreLevels();
  Logger* pLog_;
  Orbitals& orbitals_;

  // program tasks

  bool do_gw_ = false;
  bool do_bse_singlets_ = false;
  bool do_bse_triplets_ = false;
  bool do_dynamical_screening_bse_ = false;

  // options for own Vxc calculation
  std::string functional_;
  std::string grid_;

  GW::options gwopt_;
  BSE::options bseopt_;

  std::string sigma_plot_states_;
  Index sigma_plot_steps_;
  double sigma_plot_spacing_;
  std::string sigma_plot_filename_;

  // basis sets
  std::string auxbasis_name_;
  std::string dftbasis_name_;

  std::vector<QMFragment<BSE_Population> > fragments_;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GWBSE_H
