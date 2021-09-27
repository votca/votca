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

// Local VOTCA includes
#include "votca/xtp/rate_engine.h"

namespace votca {
namespace xtp {

std::ostream& operator<<(std::ostream& out, const Rate_Engine& rate_engine) {
  out << "Rate engine initialized:\n";
  out << " Ratetype:" << rate_engine.ratetype_ << "\n";
  out << " Temperature T[k] = "
      << rate_engine.temperature_ * tools::conv::hrt2ev / tools::conv::kB
      << "\n";
  Eigen::Vector3d field =
      rate_engine.field_ * tools::conv::hrt2ev / tools::conv::bohr2nm;
  out << " Electric field[V/nm](x,y,z) =" << field.x() << " " << field.y()
      << " " << field.z() << " ||F|| " << field.norm() << std::endl;
  return out;
}

Rate_Engine::PairRates Rate_Engine::Rate(const QMPair& pair,
                                         QMStateType carriertype) const {
  double charge = 0.0;
  if (carriertype == QMStateType::Electron) {
    charge = -1.0;
  } else if (carriertype == QMStateType::Hole) {
    charge = 1.0;
  }

  double reorg12 = pair.getReorg12(carriertype) + pair.getLambdaO(carriertype);
  double reorg21 = pair.getReorg21(carriertype) - pair.getLambdaO(carriertype);
  if (std::abs(reorg12) < 1e-12 || std::abs(reorg21) < 1e-12) {
    throw std::runtime_error(
        "Reorganisation energy for a pair is extremely close to zero,\n"
        " you probably forgot to import reorganisation energies into your "
        "state "
        "file.");
  }
  double dG_Field = 0.0;
  if (charge != 0.0) {
    dG_Field = charge * pair.R().dot(field_);
  }
  double dG_Site = pair.getdE12(carriertype);
  double dG = dG_Site + dG_Field;
  double J2 = pair.getJeff2(carriertype);
  PairRates result;
  if (ratetype_ == "marcus") {
    result.rate12 = Marcusrate(J2, dG, reorg12);
    result.rate21 = Marcusrate(J2, -dG, reorg21);
  } else {
    throw std::runtime_error("Only marcus rates implemented.");
  }
  return result;
}

double Rate_Engine::Marcusrate(double Jeff2, double deltaG,
                               double reorg) const {

  double hbar = tools::conv::hbar * tools::conv::ev2hrt;
  return 2 * tools::conv::Pi / hbar * Jeff2 /
         std::sqrt(4 * tools::conv::Pi * reorg * temperature_) *
         std::exp(-(deltaG - reorg) * (deltaG - reorg) /
                  (4 * reorg * temperature_));
}
}  // namespace xtp
}  // namespace votca
