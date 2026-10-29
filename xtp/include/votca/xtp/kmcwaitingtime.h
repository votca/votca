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
#ifndef VOTCA_XTP_KMCWAITINGTIME_H
#define VOTCA_XTP_KMCWAITINGTIME_H

// Standard includes
#include <memory>
#include <string>
#include <vector>

// Local VOTCA includes
#include "kmcrandom.h"
#include "kmcrates.h"

namespace votca {
namespace xtp {

/**
 * \brief How a carrier's next hop is drawn: the waiting time and the
 * destination. The pluggable part of the KMC engine (KMCEngine).
 *
 * Every carrier has its own clock. The engine asks NextTime() for the time of
 * the carrier's next candidate event; at that time it asks ChooseHop() for
 * the edge, which may return -1 (a null event, e.g. a rejected thinning
 * proposal). If the chosen destination is occupied by another carrier the
 * engine also treats the event as a null event. After every event of a
 * carrier (hop or null) the engine asks NextTime() again from the event
 * time; the clocks of the other carriers are not touched.
 *
 * This is exact if, for a carrier on site i, the candidate events form a
 * (possibly time-inhomogeneous) Poisson process with intensity c_i(t) and
 * ChooseHop returns edge e with probability w_e(t)/c_i(t), w the rates of
 * the TimeDependentRates, c_i(t) >= sum_e w_e(t). Exclusion then enters only
 * through the rejected events, which is why strategies use the full escape
 * rates of the sites and never look at the occupation.
 *
 * Contract for implementations:
 *  - Prepare() is called once, before any run, single threaded; put all
 *    precomputation there (e.g. integrated escape rates per site).
 *  - NextTime() and ChooseHop() are const and called concurrently from
 *    several runs; all randomness must come from the KMCRandom passed in, so
 *    that runs are reproducible.
 *  - NextTime() returns a time > t (or +infinity if the carrier can never
 *    move again).
 *
 * Implemented: ExponentialWaitingTime (static rates), ThinningWaitingTime
 * (any rates with upper bounds). Template for the integrated-rate inversion
 * for periodic rates (Prados, Brey, Sanchez-Rey, J. Stat. Phys. 89, 709
 * (1997); R. Paiders, BEP TU/e 2025): IntegratedRateInversion.
 */
class WaitingTimeStrategy {
 public:
  virtual ~WaitingTimeStrategy() = default;
  virtual std::string Name() const = 0;
  virtual void Prepare(const TimeDependentRates& rates) = 0;
  virtual double NextTime(Index site, double t, KMCRandom& random) const = 0;
  virtual Index ChooseHop(Index site, double t, KMCRandom& random) const = 0;
};

/**
 * \brief Time-homogeneous KMC: waiting time Exp(escape rate), destination
 * with probability w_e / escape rate. Only for static rates.
 */
class ExponentialWaitingTime final : public WaitingTimeStrategy {
 public:
  std::string Name() const override { return "exponential"; }
  void Prepare(const TimeDependentRates& rates) override;
  double NextTime(Index site, double t, KMCRandom& random) const override;
  Index ChooseHop(Index site, double t, KMCRandom& random) const override;

 private:
  const ToFGraph* graph_ = nullptr;
  std::vector<double> escape_;
  std::vector<double> cumulative_;  // per edge, cumulative within its site
};

/**
 * \brief Thinning (Lewis-Shedler) for any rates with time-independent upper
 * bounds B_e >= w_e(t): candidate events at the constant rate sum_e B_e,
 * destination proposed with probability B_e / sum B, accepted with
 * probability w_e(t)/B_e. Exact for every time dependence; the cost grows
 * with the ratio bound / mean rate.
 */
class ThinningWaitingTime final : public WaitingTimeStrategy {
 public:
  std::string Name() const override { return "thinning"; }
  void Prepare(const TimeDependentRates& rates) override;
  double NextTime(Index site, double t, KMCRandom& random) const override;
  Index ChooseHop(Index site, double t, KMCRandom& random) const override;

 private:
  const TimeDependentRates* rates_ = nullptr;
  std::vector<double> bound_escape_;
  std::vector<double> cumulative_;  // of the bounds
};

/**
 * \brief TEMPLATE, NOT IMPLEMENTED: waiting times for periodic rates by
 * inverting the integrated escape rate (BEP R. Paiders, nhKMC, Algorithms 3
 * and 4; Prados, Brey, Sanchez-Rey 1997).
 *
 * For a carrier on site i at time t0 the waiting time dt solves
 *   -ln u = int_{t0}^{t0+dt} lambda_i(s) ds,   u ~ U(0,1],
 * lambda_i(t) = TimeDependentRates::EscapeRate(i, t) > 0 with period T.
 * With Lambda_i(t) = int_0^t lambda_i, I_T = Lambda_i(T) and t0 = l T + t0*:
 *   m = floor((g + Lambda_i(t0*)) / I_T) ... (BEP eqs. 11-18)
 * so only Lambda_i on [0, T) (or [0, 2T)) and I_T are needed, precomputed per
 * site in Prepare() (quadrature + monotone PCHIP interpolant), and the
 * remainder is a bounded, monotone root search.
 *
 * ChooseHop(i, t): edge e with probability w_e(t) / lambda_i(t) at the event
 * time t (never -1, apart from rounding).
 *
 * Notes for the implementation:
 *  - Prepare(): rates.Period() must be > 0, else throw. Store per site what
 *    the inversion needs (e.g. the PCHIP coefficients of Lambda_i on [0, T]
 *    and I_T); memory is per site, not per edge.
 *  - Tolerances relative to T and I_T that double precision can reach
 *    (about 1e-14 relative, not 1e-20). The phase rem(t0, T) loses about
 *    log10(t0/T) digits for t0 >> T; document the accuracy or keep the time
 *    split into periods and phase.
 *  - Validation (see test_kmcwaitingtime.cc, disabled test cases):
 *    * static rates (Period() any > 0): identical statistics to
 *      ExponentialWaitingTime and the GRW hitting times;
 *    * uniformly scaled rates w_e(t) = s(t) w_e: hitting times are the
 *      static ones mapped through the inverse of int s (exact reference);
 *    * any periodic rates: same statistics as ThinningWaitingTime;
 *    * f -> infinity: the GRW with period-averaged rates; f -> 0: the
 *      average over the phase of the static result.
 */
class IntegratedRateInversion final : public WaitingTimeStrategy {
 public:
  struct Settings {
    Index quadrature_points = 1000;  // per period, for Lambda_i
    Index interpolation_points = 100;
    double tolerance = 1e-14;  // relative, root search
  };
  explicit IntegratedRateInversion(const Settings& settings)
      : settings_(settings) {}
  std::string Name() const override { return "integrated-rate inversion"; }
  void Prepare(const TimeDependentRates& rates) override;
  double NextTime(Index site, double t, KMCRandom& random) const override;
  Index ChooseHop(Index site, double t, KMCRandom& random) const override;

 private:
  Settings settings_;
  const TimeDependentRates* rates_ = nullptr;
};

/// strategy by name: exponential, thinning, inversion
std::unique_ptr<WaitingTimeStrategy> MakeWaitingTimeStrategy(
    const std::string& name);

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCWAITINGTIME_H
