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
  std::string key = identify();
  std::string filename =
      prop.ifExistsReturnElseThrowRuntimeError<std::string>("options_polar");
  tools::Property polar_xml;
  tools::load_property_from_xml(polar_xml, filename);
  _max_iter =
      polar_xml.ifExistsReturnElseReturnDefault(key + ".max_iter", _max_iter);
  _deltaE = polar_xml.ifExistsReturnElseReturnDefault(key + ".tolerance_energy",
                                                      _deltaE);
  _deltaD = polar_xml.ifExistsReturnElseReturnDefault(key + ".tolerance_dipole",
                                                      _deltaD);
  _exp_damp =
      polar_xml.ifExistsReturnElseReturnDefault(key + ".exp_damp", _exp_damp);
  _induce_intra_mol = polar_xml.ifExistsReturnElseReturnDefault(
      key + ".induce_intra_molecule", _induce_intra_mol);
  _openmp_threads = polar_xml.ifExistsReturnElseReturnDefault<int>(
      key + ".openmp", _openmp_threads);

  OPENMP::setMaxThreads(_openmp_threads);
  return;
}

bool PolarRegion::Converged() const {

  if (!_E_hist.filled()) {
    return false;
  }
  double Echange = _E_hist.getDiff();
  double Dchange = _D_hist.getDiff();
  std::string info = "not converged";
  bool converged = false;
  if (std::abs(Dchange) < _deltaD && std::abs(Echange) < _deltaE) {
    info = "converged";
    converged = true;
  }
  XTP_LOG_SAVE(logINFO, _log)
      << "Region:" << this->identify() << " " << this->getId() << " is " << info
      << " deltaE=" << Echange << " deltaDipole=" << Dchange << std::flush;
  return converged;
}

double PolarRegion::StaticInteraction() {
  double energy_ext = 0.0;
  for (const PolarSegment& seg : _segments) {
    for (const PolarSite site : seg) {
      energy_ext += site.Energy();
    }
  }
  double e_static = 0.0;
  eeInteractor eeinteractor;
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < i; ++j) {
      e_static += eeinteractor.InteractStatic(_segments[i], _segments[j]);
    }
  }
  if (_induce_intra_mol) {
    for (PolarSegment& seg : _segments) {
      e_static += eeinteractor.InteractStatic_IntraSegment(seg);
    }
  }
  return energy_ext + e_static;
}

void PolarRegion::PolarInteraction_scf() {
  eeInteractor eeinteractor(_exp_damp);

  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < i; ++j) {
      eeinteractor.InteractPolar_small(_segments[i], _segments[j]);
    }
  }
  if (_induce_intra_mol) {
    for (PolarSegment& seg : _segments) {
      eeinteractor.Cholesky_IntraSegment(seg);
    }
  }
}

double PolarRegion::PolarInteraction_energy() {
  double e = 0.0;
  eeInteractor eeinteractor(_exp_damp);
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < i; ++j) {
      e += eeinteractor.InteractPolar(_segments[i], _segments[j]);
    }
  }
  if (_induce_intra_mol) {
    for (PolarSegment& seg : _segments) {
      e += eeinteractor.InteractPolar_IntraSegment(seg);
    }
  }
  for (PolarSegment& seg : _segments) {
    for (PolarSite& site : seg) {
      e += site.InternalEnergy();
    }
  }

  return e;
}

void PolarRegion::CalcInducedDipoles() {
  for (PolarSegment& seg : _segments) {
    for (PolarSite& site : seg) {
      site.calc_InducedDipole();
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
        maxchange = change;
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

  ApplyInfluenceOfOtherRegions(regions);
  XTP_LOG_SAVE(logINFO, _log)
      << "Evaluating electrostatics inside region" << std::flush;
  double e_static = StaticInteraction();
  XTP_LOG_SAVE(logINFO, _log)
      << "Calculated static energy[hrt]= " << e_static << std::flush;
  XTP_LOG_SAVE(logINFO, _log)
      << "Starting SCF for classical polarisation" << std::flush;

  if (_induce_intra_mol) {
    eeInteractor eeinteractor(_exp_damp);
    for (PolarSegment& seg : _segments) {
      eeinteractor.Cholesky_IntraSegment(seg);
    }
  }
  for (int iteration = 0; iteration < _max_iter; iteration++) {
    XTP_LOG_SAVE(logINFO, _log) << identify() << " Iteration " << iteration + 1
                                << " of " << _max_iter << std::flush;
    CalcInducedDipoles();
    XTP_LOG_SAVE(logINFO, _log)
        << "Calculated induced dipoles from fields" << std::flush;
    ResetFields();
    XTP_LOG_SAVE(logINFO, _log) << "Reset old induced fields" << std::flush;
    PolarInteraction_scf();

    std::pair<bool, double> d_converged = DipolesConverged();
    bool e_converged = false;

    bool converged = d_converged.first;
    if (converged || (iteration == (_max_iter - 1))) {
      if (converged) {
        XTP_LOG_SAVE(logINFO, _log)
            << "SCF calculation converged after " << iteration + 1
            << " iterations." << std::flush;
      } else {
        XTP_LOG_SAVE(logINFO, _log)
            << "WARNING: SCF calculation not converged after " << iteration + 1
            << " iterations!!!" << std::flush;
      }
      _D_hist.push_back(d_converged.second);
      break;
    }
  }
  ResetFields();
  double e_polar = PolarInteraction_energy();
  XTP_LOG_SAVE(logINFO, _log)
      << " E_polar[hrt]= " << e_polar << " E_static[hrt]= " << e_static
      << " E_total[hrt]= " << e_static + e_polar << std::flush;
  _E_hist.push_back(e_static + e_polar);
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
