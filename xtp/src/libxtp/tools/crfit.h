/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_CRFIT_H
#define VOTCA_XTP_CRFIT_H

// Local VOTCA includes
#include "votca/xtp/chargeresponsefit.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

/**
 * \brief Charge-response parametrization of a molecule with xtpdft: the
 * charge-response matrix Pi (.crm) and the residual polarizabilities (.mps,
 * from a template with the permanent multipoles), from central finite
 * differences of SCFs in Hirshfeld-localized potentials and uniform fields
 * (see ChargeResponseFit).
 */
class CRFit final : public QMTool {
 public:
  std::string Identify() const final { return "crfit"; }

 protected:
  void ParseOptions(const tools::Property& options) final;
  bool Run() final;

 private:
  void WriteCrm(const QMMolecule& mol,
                const ChargeResponseFit::Result& r) const;
  void WriteJson(const QMMolecule& mol,
                 const ChargeResponseFit::Result& r) const;
  void WriteMps(const QMMolecule& mol,
                const ChargeResponseFit::Result& r) const;

  mutable Logger log_;
  tools::Property dftoptions_;
  std::string xyzfile_;
  std::string mpsfile_;
  std::string output_;
  double step_ = 1e-3;
  double scf_error_ = 1e-9;
  bool inversion_ = true;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CRFIT_H
