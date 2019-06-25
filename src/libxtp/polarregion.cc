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

#include "votca/xtp/eeinteractor.h"
#include <votca/xtp/dipoledipoleinteraction.h>
#include <votca/xtp/polarregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/staticregion.h>

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
  _deltaD = polar_xml.ifExistsReturnElseReturnDefault(key + ".tolerance_dipole",
                                                      _deltaD);
  _deltaE = polar_xml.ifExistsReturnElseReturnDefault(key + ".tolerance_energy",
                                                      _deltaE);
  _exp_damp =
      polar_xml.ifExistsReturnElseReturnDefault(key + ".exp_damp", _exp_damp);
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
  std::string info = "not converged";
  bool converged = false;
  if (std::abs(Echange) < _deltaE) {
    info = "converged";
    converged = true;
  }
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << "Region:" << this->identify() << " " << this->getId()
      << " is " << info << " deltaE=" << Echange << std::flush;
  return converged;
}

void PolarRegion::StaticInteraction() {

  eeInteractor eeinteractor;
#pragma omp parallel for
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < size(); ++j) {
      if (i == j) {
        continue;
      }
      eeinteractor.ApplyStaticField(_segments[i], _segments[j]);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < size(); ++i) {
    eeinteractor.ApplyStaticField_IntraSegment(_segments[i]);
  }
}

double PolarRegion::Energy() const {
  double e = 0.0;
  eeInteractor eeinteractor(_exp_damp);
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < i; ++j) {
      e += eeinteractor.CalcStaticEnergy(_segments[i], _segments[j]);
      e += eeinteractor.CalcPolarEnergy(_segments[i], _segments[j]);
    }
  }

  for (const PolarSegment& seg : _segments) {
    e += eeinteractor.CalcPolarEnergy_IntraSegment(seg);
    e += eeinteractor.CalcStaticEnergy_IntraSegment(seg);
  }
  for (const PolarSegment& seg : _segments) {
    for (const PolarSite& site : seg) {
      e += site.InternalEnergy();
    }
  }
  return e;
}

void PolarRegion::Reset() {
  for (PolarSegment& seg : _segments) {
    for (PolarSite& site : seg) {
      site.Reset();
    }
  }
}

void PolarRegion::Evaluate(std::vector<std::unique_ptr<Region> >& regions) {
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Evaluating:" << this->identify() << " "
      << this->getId() << std::flush;

  ApplyInfluenceOfOtherRegions(regions);
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Evaluating electrostatics inside region"
      << std::flush;
  StaticInteraction();
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Calculated static interaction in region"
      << std::flush;

  int dof_polarisation = 0;
  for (const PolarSegment& seg : _segments) {
    dof_polarisation += seg.size() * 3;
  }
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Starting Solving for classical polarisation with "
      << dof_polarisation << " degrees of freedom." << std::flush;

  Eigen::VectorXd initial_induced_dipoles =
      Eigen::VectorXd::Zero(dof_polarisation);

  if (!_E_hist.filled()) {
    eeInteractor interactor(_exp_damp);
    int index = 0;
    for (PolarSegment& seg : _segments) {
      initial_induced_dipoles.segment(index, 3 * seg.size()) =
          interactor.Cholesky_IntraSegment(seg);
      index += 3 * seg.size();
    }
  } else {
    int index = 0;
    for (PolarSegment& seg : _segments) {
      for (const PolarSite& site : seg) {
        initial_induced_dipoles.segment<3>(index) = site.Induced_Dipole();
        index += 3;
      }
    }
  }
  Eigen::VectorXd b = Eigen::VectorXd::Zero(dof_polarisation);
  int index = 0;
  for (PolarSegment& seg : _segments) {
    for (const PolarSite& site : seg) {
      b.segment<3>(index) = -site.V().segment<3>(1);
      index += 3;
    }
  }

  eeInteractor interactor(_exp_damp);
  DipoleDipoleInteraction A(interactor, _segments);
  Eigen::ConjugateGradient<DipoleDipoleInteraction, Eigen::Lower | Eigen::Upper>
      cg;
  cg.setMaxIterations(_max_iter);
  cg.setTolerance(_deltaD);
  cg.compute(A);
  Eigen::VectorXd x = cg.solveWithGuess(b, initial_induced_dipoles);

  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " CG: #iterations: " << cg.iterations()
      << ", estimated error: " << cg.error() << std::endl;
  index = 0;
  for (PolarSegment& seg : _segments) {
    for (PolarSite& site : seg) {
      site.setInduced_Dipole(x.segment<3>(index));
      index += 3;
    }
  }

  double e_total = Energy();
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " E_total[hrt]= " << e_total << std::flush;
  _E_hist.push_back(e_total);
  return;
}

void PolarRegion::InteractwithQMRegion(const QMRegion& region) {
  region.ApplyQMFieldToPolarSegments(_segments);
}
void PolarRegion::InteractwithPolarRegion(const PolarRegion& region) {
#pragma omp parallel for
  for (unsigned i = 0; i < _segments.size(); i++) {
    PolarSegment& pseg1 = _segments[i];
    eeInteractor ee;
    for (const PolarSegment& pseg2 : region) {
      ee.ApplyStaticField(pseg2, pseg1);
      ee.ApplyInducedField(pseg2, pseg1);
    }
  }
}

void PolarRegion::InteractwithStaticRegion(const StaticRegion& region) {
#pragma omp parallel for
  for (unsigned i = 0; i < _segments.size(); i++) {
    PolarSegment& pseg = _segments[i];
    eeInteractor ee;
    for (const StaticSegment& sseg : region) {
      ee.ApplyStaticField(sseg, pseg);
    }
  }

  return;
}

}  // namespace xtp
}  // namespace votca
