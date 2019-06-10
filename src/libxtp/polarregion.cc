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

#include <votca/xtp/polarregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/staticregion.h>

#include "votca/xtp/eeinteractor.h"

namespace votca {
namespace xtp {

void PolarRegion::Initialize(const tools::Property& prop) {

  std::string filename =
      prop.ifExistsReturnElseThrowRuntimeError<std::string>("options");
  tools::Property polar_xml;
  tools::load_property_from_xml(polar_xml, filename);
  _max_iter = polar_xml.ifExistsReturnElseReturnDefault("max_iter", _max_iter);
  _deltaE =
      polar_xml.ifExistsReturnElseReturnDefault("tolerance_energy", _deltaE);
  _deltaD =
      polar_xml.ifExistsReturnElseReturnDefault("tolerance_dipole", _deltaD);
  _exp_damp = polar_xml.ifExistsReturnElseReturnDefault("exp_damp", _exp_damp);
  _induce_intra_mol = polar_xml.ifExistsReturnElseReturnDefault(
      "induce_intra_molecule", _induce_intra_mol);

  return;
}

bool PolarRegion::Converged() const {

  double Echange = std::abs(_E_hist.getDiff());
  double Dchange = std::abs(_D_hist.getDiff());
  std::string info = "not converged";
  bool converged = false;
  if (Dchange < _deltaD && Echange < _deltaE) {
    info = "converged";
    converged = true;
  }
  XTP_LOG_SAVE(logINFO, _log)
      << "Region:" << this->identify() << " " << this->getId() << " is " << info
      << " deltaE=" << Echange << " deltaDipole=" << Dchange << std::flush;
  return converged;
}

double PolarRegion::StaticInteraction() {
  double static_energy = 0.0;
  eeInteractor eeinteractor;
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < i; ++j) {
      static_energy += eeinteractor.InteractStatic(_segments[i], _segments[j]);
    }
  }
  if (_induce_intra_mol) {
    for (PolarSegment& seg : _segments) {
      static_energy += eeinteractor.InteractStatic_IntraSegment(seg);
    }
  }
  return static_energy;
}

double PolarRegion::PolarInteraction() {
  double polar_energy = 0.0;
  eeInteractor eeinteractor(_exp_damp);
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < i; ++j) {
      polar_energy += eeinteractor.InteractPolar(_segments[i], _segments[j]);
    }
  }
  if (_induce_intra_mol) {
    for (PolarSegment& seg : _segments) {
      polar_energy += eeinteractor.InteractPolar_IntraSegment(seg);
    }
  }
  return polar_energy;
}

void PolarRegion::CalcInducedDipoles() {
  for (PolarSegment& seg : _segments) {
    for (PolarSite& site : seg) {
      site.calcDIIS_InducedDipole();
    }
  }
}

void PolarRegion::ResetFields() {
  for (PolarSegment& seg : _segments) {
    for (PolarSite& site : seg) {
      site.ResetInduction();
    }
  }
}

std::pair<bool, double> PolarRegion::DipolesConverged() const {
  int converged_sites = 0;
  int sites = 0;
  double maxchange = 0.0;
  int segid = 0;
  int siteid = 0;
  for (const PolarSegment& seg : _segments) {
    sites += seg.size();
    for (const PolarSite& site : seg) {
      double change = site.DipoleChange();
      if (change > maxchange) {
        change = maxchange;
        segid = seg.getId();
        siteid = site.getId();
      }
      if (change < _deltaD) {
        converged_sites++;
      }
    }
  }
  std::pair<bool, double> result;
  result.first = (converged_sites == sites);
  result.second = maxchange;
  double percent = 100 * double(converged_sites) / double(sites);
  XTP_LOG_SAVE(logINFO, _log)
      << percent << "% of sites converged. Max dipole change:" << maxchange
      << " segment:" << segid << " site:" << siteid << std::flush;
  return result;
}

void PolarRegion::Evaluate(std::vector<std::unique_ptr<Region> >& regions) {
  XTP_LOG_SAVE(logINFO, _log) << "Evaluating:" << this->identify() << " "
                              << this->getId() << std::flush;
  ResetRegion();
  XTP_LOG_SAVE(logINFO, _log)
      << "Removed all previous values from region" << std::flush;
  ApplyInfluenceOfOtherRegions(regions);
  XTP_LOG_SAVE(logINFO, _log)
      << "Evaluating electrostatics inside region" << std::flush;
  double energy = StaticInteraction();
  XTP_LOG_SAVE(logINFO, _log)
      << "Starting SCF for classical polarisation" << std::flush;

  double e_old = 0.0;
  for (int iteration = 0; iteration < _max_iter; iteration++) {
    XTP_LOG_SAVE(logINFO, _log)
        << "Iteration " << iteration + 1 << " of " << _max_iter << std::flush;
    CalcInducedDipoles();
    XTP_LOG_SAVE(logINFO, _log)
        << "Calculated induced dipoles from fields" << std::flush;
    ResetFields();
    XTP_LOG_SAVE(logINFO, _log) << "Reset old induced fields" << std::flush;
    double polar_energy = PolarInteraction();
    XTP_LOG_SAVE(logINFO, _log)
        << "Calculated polar interactions" << std::flush;

    std::pair<bool, double> d_converged = DipolesConverged();
    bool e_converged = false;
    double deltaE = std::abs(polar_energy - e_old);
    e_old = polar_energy;
    if (iteration > 0) {
      XTP_LOG_SAVE(logINFO, _log)
          << "Change from last iteration DeltaE [Ha]:" << deltaE << std::flush;
      if (deltaE < _deltaE) {
        e_converged = true;
        XTP_LOG_SAVE(logINFO, _log)
            << "Energy converged to" << _deltaE << std::flush;
      }
    }
    bool converged = d_converged.first && e_converged;
    if (converged || (iteration == (_max_iter - 1))) {
      energy += polar_energy;
      _E_hist.push_back(energy);
      _D_hist.push_back(d_converged.second);
      if (converged) {
        XTP_LOG_SAVE(logINFO, _log)
            << "SCF calculation converged after " << iteration + 1
            << " iterations." << std::flush;
      } else {
        XTP_LOG_SAVE(logINFO, _log)
            << "WARNING: SCF calculation not converged after " << iteration + 1
            << " iterations!!!" << std::flush;
      }

      break;
    }
  }

  return;
}

void PolarRegion::InteractwithQMRegion(const QMRegion& region) {
  region.ApplyQMFieldToClassicSegments(_segments);
}
void PolarRegion::InteractwithPolarRegion(const PolarRegion& region) {
#pragma omp parallel for
  for (int i = 0; i < int(_segments.size()); i++) {
    PolarSegment& pseg1 = _segments[i];
    eeInteractor ee;

    for (const PolarSegment& pseg2 : region) {
      ee.InteractStatic(pseg2, pseg1);
      ee.InteractPolar(pseg2, pseg1);
    }
  }
}
void PolarRegion::InteractwithStaticRegion(const StaticRegion& region) {
#pragma omp parallel for
  for (int i = 0; i < int(_segments.size()); i++) {
    PolarSegment& pseg = _segments[i];
    eeInteractor ee;

    for (const StaticSegment& sseg : region) {
      ee.InteractStatic(sseg, pseg);
    }
  }

  return;
}

}  // namespace xtp
}  // namespace votca
