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
#include <numeric>
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

  return;
}

bool PolarRegion::Converged() const {

  if (!_E_hist.filled()) {
    return false;
  }
  double Echange = _E_hist.getDiff().Etotal();
  std::string info = "not converged";
  bool converged = false;
  if (std::abs(Echange) < _deltaE) {
    info = "converged";
    converged = true;
  }
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Region:" << this->identify() << " " << this->getId()
      << " is " << info << " deltaE=" << Echange << std::flush;
  return converged;
}

double PolarRegion::StaticInteraction() {

  eeInteractor eeinteractor;
  double e = 0.0;
#pragma omp parallel for reduction(+ : e)
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < size(); ++j) {
      if (i == j) {
        continue;
      }
      PolarSegment& seg1 = _segments[i];
      const PolarSegment& seg2 = _segments[j];
      e += eeinteractor.ApplyStaticField<PolarSegment, Estatic::noE_V>(seg2,
                                                                       seg1);
    }
  }

  return e / 2.0;
}

eeInteractor::E_terms PolarRegion::PolarEnergy() const {
#pragma omp declare reduction(CustomPlus              \
                              : eeInteractor::E_terms \
                              : omp_out += omp_in)

  eeInteractor eeinteractor(_exp_damp);

  eeInteractor::E_terms terms;

#pragma omp parallel for reduction(CustomPlus : terms)
  for (int i = 0; i < size(); ++i) {
    for (int j = 0; j < i; ++j) {
      terms += eeinteractor.CalcPolarEnergy(_segments[i], _segments[j]);
    }
  }

#pragma omp parallel for reduction(CustomPlus : terms)
  for (int i = 0; i < size(); ++i) {
    terms.E_indu_indu() +=
        eeinteractor.CalcPolarEnergy_IntraSegment(_segments[i]);
  }

#pragma omp parallel for reduction(CustomPlus : terms)
  for (int i = 0; i < size(); ++i) {
    for (const PolarSite& site : _segments[i]) {
      terms.E_internal() += site.InternalEnergy();
    }
  }
  return terms;
}

double PolarRegion::PolarEnergy_extern() const {
  double e = 0.0;
#pragma omp parallel for reduction(+ : e)
  for (int i = 0; i < size(); ++i) {
    for (const PolarSite& site : _segments[i]) {
      e += site.deltaQ_V_ext();
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

void PolarRegion::AppendResult(tools::Property& prop) const {
  const Energy_terms& e = this->_E_hist.back();
  prop.add("E_static", std::to_string(e.Estatic() * tools::conv::hrt2ev));
  prop.add("E_polar", std::to_string(e.Epolar() * tools::conv::hrt2ev));
  prop.add("E_total", std::to_string(e.Etotal() * tools::conv::hrt2ev));
}

void PolarRegion::Evaluate(std::vector<std::unique_ptr<Region> >& regions) {

  std::vector<double> energies = ApplyInfluenceOfOtherRegions(regions);
  Energy_terms e_contrib;
  e_contrib.E_static_ext() =
      std::accumulate(energies.begin(), energies.end(), 0.0);
  XTP_LOG_SAVE(logINFO, _log) << TimeStamp()
                              << " Calculated static-static and polar-static "
                                 "interaction with other regions"
                              << std::flush;
  e_contrib.E_static_static() = StaticInteraction();
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Calculated static interaction in region "
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
      const Eigen::Vector3d V = site.V() + site.V_noE();
      b.segment<3>(index) = -V;
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
      << ", estimated error: " << cg.error() << std::flush;

  if (cg.info() == Eigen::ComputationInfo::NoConvergence) {
    _info = false;
    _errormsg = "PCG iterations did not converge";
    return;
  }

  index = 0;
  for (PolarSegment& seg : _segments) {
    for (PolarSite& site : seg) {
      site.setInduced_Dipole(x.segment<3>(index));
      index += 3;
    }
  }

  e_contrib.addInternalPolarContrib(PolarEnergy());
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Calculated polar interaction in region" << std::flush;
  e_contrib.E_polar_ext() = PolarEnergy_extern();
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Calculated polar interaction with other regions"
      << std::flush;

  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "   Internal static energy [hrt]= " << e_contrib.E_static_static()
      << std::flush;
  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "   External static energy [hrt]= " << e_contrib.E_static_ext()
      << std::flush;
  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "  Total static energy [hrt]= " << e_contrib.Estatic() << std::flush;

  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "   internal dQ-dQ energy [hrt]= " << e_contrib.E_indu_indu()
      << std::flush;

  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "   internal Q-dQ energy [hrt]= " << e_contrib.E_indu_stat()
      << std::flush;

  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "   Internal energy [hrt]= " << e_contrib.E_internal() << std::flush;

  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "   External polar energy [hrt]= " << e_contrib.E_polar_ext()
      << std::flush;

  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10)
      << "  Total polar energy [hrt]= " << e_contrib.Epolar() << std::flush;

  XTP_LOG_SAVE(logINFO, _log)
      << std::setprecision(10) << " Total energy [hrt]= " << e_contrib.Etotal()
      << std::flush;
  _E_hist.push_back(e_contrib);
  return;
}

double PolarRegion::InteractwithQMRegion(const QMRegion& region) {
  // QMregions always have lower ids than other regions
  region.ApplyQMFieldToPolarSegments(_segments);
  return 0.0;
}
double PolarRegion::InteractwithPolarRegion(const PolarRegion& region) {
  bool noE_V = true;
  if (this->getId() < region.getId()) {
    noE_V = false;
  }

  double e = 0;
#pragma omp parallel for reduction(+ : e)
  for (unsigned i = 0; i < _segments.size(); i++) {
    PolarSegment& pseg1 = _segments[i];
    double e_thread = 0.0;
    eeInteractor ee;
    for (const PolarSegment& pseg2 : region) {
      if (noE_V) {
        ee.ApplyStaticField<PolarSegment, Estatic::noE_V>(pseg2, pseg1);
        ee.ApplyInducedField<Estatic::noE_V>(pseg2, pseg1);
      } else {
        e_thread += ee.ApplyStaticField<PolarSegment, Estatic::V>(pseg2, pseg1);
        e_thread += ee.ApplyInducedField<Estatic::V>(pseg2, pseg1);
      }
      e += e_thread;
    }
  }
  return e;
}

double PolarRegion::InteractwithStaticRegion(const StaticRegion& region) {
  // Static regions always have higher ids than other regions

  double e = 0.0;
#pragma omp parallel for reduction(+ : e)
  for (unsigned i = 0; i < _segments.size(); i++) {
    double e_thread = 0.0;
    PolarSegment& pseg = _segments[i];
    eeInteractor ee;
    for (const StaticSegment& sseg : region) {
      e_thread += ee.ApplyStaticField<StaticSegment, Estatic::V>(sseg, pseg);
    }
    e += e_thread;
  }

  return e;
}

void PolarRegion::WriteToCpt(CheckpointWriter& w) const {
  MMRegion<PolarSegment>::WriteToCpt(w);
}

void PolarRegion::ReadFromCpt(CheckpointReader& r) {
  MMRegion<PolarSegment>::ReadFromCpt(r);
}

}  // namespace xtp
}  // namespace votca
