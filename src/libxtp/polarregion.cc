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

// Standard includes
#include <iomanip>
#include <numeric>

// Local VOTCA includes
#include "votca/xtp/dipoledipoleinteraction.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/qmregion.h"
#include "votca/xtp/staticregion.h"

namespace votca {
namespace xtp {

void PolarRegion::Initialize(const tools::Property& prop) {
<<<<<<< HEAD
  max_iter_ = prop.get( "polar.max_iter").as<Index>();
  deltaD_ = prop.get("polar.tolerance_dipole").as<double>();
  deltaE_ =prop.get("polar.tolerance_energy").as<double>();
  exp_damp_ = prop.get("polar.exp_damp").as<double>();
=======
  std::string key = identify();
  tools::Property options_polar;
  options_polar.add("polar", "");
  tools::Property& prop_polar = options_polar.get("polar");
  prop_polar = prop;

  max_iter_ = options_polar.ifExistsReturnElseReturnDefault(key + ".max_iter",
                                                            max_iter_);
  deltaD_ = options_polar.ifExistsReturnElseReturnDefault(
      key + ".tolerance_dipole", deltaD_);
  deltaE_ = options_polar.ifExistsReturnElseReturnDefault(
      key + ".tolerance_energy", deltaE_);
  exp_damp_ = options_polar.ifExistsReturnElseReturnDefault(key + ".exp_damp",
                                                            exp_damp_);
>>>>>>> 378a7d7bf2d7516765d2d830b610a907481a15ac
}

bool PolarRegion::Converged() const {

  if (!E_hist_.filled()) {
    return false;
  }
  double Echange = E_hist_.getDiff().Etotal();
  std::string info = "not converged";
  bool converged = false;
  if (std::abs(Echange) < deltaE_) {
    info = "converged";
    converged = true;
  }
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Region:" << this->identify() << " " << this->getId()
      << " is " << info << " deltaE=" << Echange << std::flush;
  return converged;
}

double PolarRegion::StaticInteraction() {

  eeInteractor eeinteractor;
  double e = 0.0;
#pragma omp parallel for reduction(+ : e)
  for (Index i = 0; i < size(); ++i) {
    for (Index j = 0; j < size(); ++j) {
      if (i == j) {
        continue;
      }
      PolarSegment& seg1 = segments_[i];
      const PolarSegment& seg2 = segments_[j];
      e += eeinteractor.ApplyStaticField<PolarSegment, Estatic::noE_V>(seg2,
                                                                       seg1);
    }
  }

  return 0.5 * e;
}

eeInteractor::E_terms PolarRegion::PolarEnergy() const {
#pragma omp declare reduction(CustomPlus              \
                              : eeInteractor::E_terms \
                              : omp_out += omp_in)

  eeInteractor eeinteractor(exp_damp_);

  eeInteractor::E_terms terms;

#pragma omp parallel for reduction(CustomPlus : terms)
  for (Index i = 0; i < size(); ++i) {
    for (Index j = 0; j < i; ++j) {
      terms += eeinteractor.CalcPolarEnergy(segments_[i], segments_[j]);
    }
  }

#pragma omp parallel for reduction(CustomPlus : terms)
  for (Index i = 0; i < size(); ++i) {
    terms.E_indu_indu() +=
        eeinteractor.CalcPolarEnergy_IntraSegment(segments_[i]);
  }

#pragma omp parallel for reduction(CustomPlus : terms)
  for (Index i = 0; i < size(); ++i) {
    for (const PolarSite& site : segments_[i]) {
      terms.E_internal() += site.InternalEnergy();
    }
  }
  return terms;
}

double PolarRegion::PolarEnergy_extern() const {
  double e = 0.0;
#pragma omp parallel for reduction(+ : e)
  for (Index i = 0; i < size(); ++i) {
    for (const PolarSite& site : segments_[i]) {
      e += site.deltaQ_V_ext();
    }
  }
  return e;
}

void PolarRegion::Reset() {
  for (PolarSegment& seg : segments_) {
    for (PolarSite& site : seg) {
      site.Reset();
    }
  }
}

void PolarRegion::AppendResult(tools::Property& prop) const {
  const Energy_terms& e = this->E_hist_.back();
  prop.add("E_static", std::to_string(e.Estatic() * tools::conv::hrt2ev));
  prop.add("E_polar", std::to_string(e.Epolar() * tools::conv::hrt2ev));
  prop.add("E_total", std::to_string(e.Etotal() * tools::conv::hrt2ev));
}

Index PolarRegion::CalcPolDoF() const {
  Index dof_polarization = 0;
  for (const PolarSegment& seg : segments_) {
    dof_polarization += seg.size() * 3;
  }
  return dof_polarization;
}

Eigen::VectorXd PolarRegion::CalcInducedDipoleInsideSegments() const {
  Index dof_polarization = CalcPolDoF();
  Eigen::VectorXd initial_induced_dipoles =
      Eigen::VectorXd::Zero(dof_polarization);
  eeInteractor interactor(exp_damp_);
  Index index = 0;
  for (const PolarSegment& seg : segments_) {
    initial_induced_dipoles.segment(index, 3 * seg.size()) =
        interactor.Cholesky_IntraSegment(seg);
    index += 3 * seg.size();
  }
  return initial_induced_dipoles;
}

Eigen::VectorXd PolarRegion::ReadInducedDipolesFromLastIteration() const {
  Index dof_polarization = CalcPolDoF();
  Eigen::VectorXd last_induced_dipoles =
      Eigen::VectorXd::Zero(dof_polarization);
  Index index = 0;
  for (const PolarSegment& seg : segments_) {
    for (const PolarSite& site : seg) {
      last_induced_dipoles.segment<3>(index) = site.Induced_Dipole();
      index += 3;
    }
  }
  return last_induced_dipoles;
}

void PolarRegion::WriteInducedDipolesToSegments(const Eigen::VectorXd& x) {
  Index index = 0;
  for (PolarSegment& seg : segments_) {
    for (PolarSite& site : seg) {
      site.setInduced_Dipole(x.segment<3>(index));
      index += 3;
    }
  }
}

Eigen::VectorXd PolarRegion::CalcInducedDipolesViaPCG(
    const Eigen::VectorXd& initial_guess) {
  Eigen::VectorXd b = Eigen::VectorXd::Zero(initial_guess.size());
  Index index = 0;
  for (PolarSegment& seg : segments_) {
    for (const PolarSite& site : seg) {
      auto V = site.V() + site.V_noE();
      b.segment<3>(index) = -V;
      index += 3;
    }
  }
  eeInteractor interactor(exp_damp_);
  DipoleDipoleInteraction A(interactor, segments_);
  Eigen::ConjugateGradient<DipoleDipoleInteraction, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg;
  cg.setMaxIterations(max_iter_);
  cg.setTolerance(deltaD_);
  cg.compute(A);
  Eigen::VectorXd x = cg.solveWithGuess(b, initial_guess);

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " CG: #iterations: " << cg.iterations()
      << ", estimated error: " << cg.error() << std::flush;

  if (cg.info() == Eigen::ComputationInfo::NoConvergence) {
    info_ = false;
    errormsg_ = "PCG iterations did not converge";
  }

  if (cg.info() == Eigen::ComputationInfo::NumericalIssue) {
    info_ = false;
    errormsg_ = "PCG had a numerical issue";
  }
  return x;
}

void PolarRegion::Evaluate(std::vector<std::unique_ptr<Region>>& regions) {

  std::vector<double> energies = ApplyInfluenceOfOtherRegions(regions);
  Energy_terms e_contrib;
  e_contrib.E_static_ext() =
      std::accumulate(energies.begin(), energies.end(), 0.0);
  XTP_LOG(Log::info, log_) << TimeStamp()
                           << " Calculated static-static and polar-static "
                              "interaction with other regions"
                           << std::flush;
  e_contrib.E_static_static() = StaticInteraction();
  XTP_LOG(Log::info, log_) << TimeStamp()
                           << " Calculated static interaction in region "
                           << std::flush;
  Index dof_polarization = CalcPolDoF();
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Starting Solving for classical polarization with "
      << dof_polarization << " degrees of freedom." << std::flush;

  Eigen::VectorXd initial_induced_dipoles;
  if (!E_hist_.filled() || segments_.size() == 1) {
    initial_induced_dipoles = CalcInducedDipoleInsideSegments();
  } else {
    initial_induced_dipoles = ReadInducedDipolesFromLastIteration();
  }

  Eigen::VectorXd x;  // if only one segment
  // it is solved exactly through the initial guess
  if (segments_.size() != 1) {
    x = CalcInducedDipolesViaPCG(initial_induced_dipoles);
    if (!info_) {
      return;
    }
  } else {
    x = initial_induced_dipoles;
  }

  WriteInducedDipolesToSegments(x);

  e_contrib.addInternalPolarContrib(PolarEnergy());
  XTP_LOG(Log::info, log_) << TimeStamp()
                           << " Calculated polar interaction in region"
                           << std::flush;
  e_contrib.E_polar_ext() = PolarEnergy_extern();
  XTP_LOG(Log::info, log_) << TimeStamp()
                           << " Calculated polar interaction with other regions"
                           << std::flush;

  XTP_LOG(Log::info, log_) << std::setprecision(10)
                           << "   Internal static energy [hrt]= "
                           << e_contrib.E_static_static() << std::flush;
  XTP_LOG(Log::info, log_) << std::setprecision(10)
                           << "   External static energy [hrt]= "
                           << e_contrib.E_static_ext() << std::flush;
  XTP_LOG(Log::error, log_)
      << std::setprecision(10)
      << "  Total static energy [hrt]= " << e_contrib.Estatic() << std::flush;

  XTP_LOG(Log::info, log_) << std::setprecision(10)
                           << "   internal dQ-dQ energy [hrt]= "
                           << e_contrib.E_indu_indu() << std::flush;

  XTP_LOG(Log::info, log_) << std::setprecision(10)
                           << "   internal Q-dQ energy [hrt]= "
                           << e_contrib.E_indu_stat() << std::flush;

  XTP_LOG(Log::info, log_) << std::setprecision(10)
                           << "   Internal energy [hrt]= "
                           << e_contrib.E_internal() << std::flush;

  XTP_LOG(Log::info, log_) << std::setprecision(10)
                           << "   External polar energy [hrt]= "
                           << e_contrib.E_polar_ext() << std::flush;

  XTP_LOG(Log::error, log_)
      << std::setprecision(10)
      << "  Total polar energy [hrt]= " << e_contrib.Epolar() << std::flush;

  XTP_LOG(Log::error, log_)
      << std::setprecision(10) << " Total energy [hrt]= " << e_contrib.Etotal()
      << std::flush;
  E_hist_.push_back(e_contrib);
  return;
}

double PolarRegion::InteractwithQMRegion(const QMRegion& region) {
  // QMregions always have lower ids than other regions
  region.ApplyQMFieldToPolarSegments(segments_);
  return 0.0;
}
double PolarRegion::InteractwithPolarRegion(const PolarRegion& region) {
  bool noE_V = true;
  if (this->getId() < region.getId()) {
    noE_V = false;
  }

  double e = 0;
#pragma omp parallel for reduction(+ : e)
  for (Index i = 0; i < Index(segments_.size()); i++) {
    PolarSegment& pseg1 = segments_[i];
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
  for (Index i = 0; i < Index(segments_.size()); i++) {
    double e_thread = 0.0;
    PolarSegment& pseg = segments_[i];
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
