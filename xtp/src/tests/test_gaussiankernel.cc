/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE gaussiankernel_test

// Standard includes
#include <cmath>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/gaussiankernel.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(gaussiankernel_test)

BOOST_AUTO_TEST_CASE(kernel_functions) {
  const double S = 1.3;
  for (double r : {0.5, 1.0, 2.0, 4.0, 9.0}) {
    GaussianKernel k(r, S);
    BOOST_CHECK_CLOSE(k.B0, std::erf(r / S) / r, 1e-12);
    // B1 = -(1/r) dB0/dr, B2 = -(1/r) dB1/dr by central differences
    const double h = 1e-5;
    GaussianKernel kp(r + h, S), km(r - h, S);
    BOOST_CHECK_CLOSE(k.B1, -(kp.B0 - km.B0) / (2 * h) / r, 1e-6);
    BOOST_CHECK_CLOSE(k.B2, -(kp.B1 - km.B1) / (2 * h) / r, 1e-6);
  }
  // point limit
  GaussianKernel p(2.0, 0.0);
  BOOST_CHECK_CLOSE(p.B0, 0.5, 1e-12);
  BOOST_CHECK_CLOSE(p.B1, 0.125, 1e-12);
  BOOST_CHECK_CLOSE(p.B2, 3.0 / 32.0, 1e-12);
  // far from the Gaussian it is the point kernel
  GaussianKernel f(20.0, 1.0);
  BOOST_CHECK_CLOSE(f.B2, 3.0 / std::pow(20.0, 5), 1e-10);
  // the short-range series joins the closed form
  for (double x : {0.05 * (1 - 1e-7), 0.05 * (1 + 1e-7)}) {
    GaussianKernel a(0.05 * S, S);
    GaussianKernel b(x * S, S);
    BOOST_CHECK_CLOSE(a.B0, b.B0, 1e-5);
    BOOST_CHECK_CLOSE(a.B1, b.B1, 1e-5);
    BOOST_CHECK_CLOSE(a.B2, b.B2, 1e-5);
  }
  GaussianKernel z(0.0, S);
  BOOST_CHECK_CLOSE(z.B0, 2.0 / std::sqrt(M_PI) / S, 1e-12);
}

namespace {
PolarSegment OneSite(Index id, const Eigen::Vector3d& pos, double alpha,
                     bool cr, double width) {
  PolarSegment seg("seg" + std::to_string(id), id);
  PolarSite site(0, "C", pos);
  site.setpolarization(alpha * Eigen::Matrix3d::Identity());
  if (cr) {
    site.setChargeResponse(width, false);
  }
  seg.push_back(site);
  return seg;
}
}  // namespace

BOOST_AUTO_TEST_CASE(induced_fields_and_energy) {
  // site 1 (charge-response) with an induced charge and dipole acting on
  // site 2 (charge-response or Thole)
  for (bool cr2 : {true, false}) {
    const double s1 = 1.2, s2 = 0.9;
    PolarSegment seg1 =
        OneSite(1, Eigen::Vector3d(0.1, -0.3, 0.2), 9.0, true, s1);
    seg1[0].setInduced_Charge(0.37);
    seg1[0].setInduced_Dipole(Eigen::Vector3d(0.2, -0.5, 0.8));
    const Eigen::Vector3d R2(2.3, 1.1, -1.7);
    auto at = [&](const Eigen::Vector3d& pos) {
      PolarSegment seg2 = OneSite(2, pos, 6.0, cr2, s2);
      eeInteractor ee;
      ee.ApplyInducedField<Estatic::noE_V>(seg1, seg2);
      return seg2;
    };
    PolarSegment seg2 = at(R2);
    const Eigen::Vector3d r = R2 - seg1[0].getPos();
    const double S = cr2 ? std::hypot(s1, s2) : s1;
    GaussianKernel k(r.norm(), S);
    if (cr2) {
      const double phi_ref =
          0.37 * k.B0 + k.B1 * seg1[0].Induced_Dipole().dot(r);
      BOOST_CHECK_CLOSE(seg2[0].Phi_noE(), phi_ref, 1e-10);
      // the field accumulated is the gradient of that potential
      const double h = 1e-5;
      for (Index c = 0; c < 3; ++c) {
        Eigen::Vector3d d = Eigen::Vector3d::Zero();
        d(c) = h;
        const double fd =
            (at(R2 + d)[0].Phi_noE() - at(R2 - d)[0].Phi_noE()) / (2 * h);
        BOOST_CHECK_CLOSE(seg2[0].V_noE()(c), fd, 1e-5);
      }
    } else {
      BOOST_CHECK_EQUAL(seg2[0].Phi_noE(), 0.0);
    }
    // pair energy = induced charge and dipole of site 2 in the potential
    // and field of site 1's induced charge and dipole
    seg2[0].setInduced_Charge(cr2 ? -0.21 : 0.0);
    seg2[0].setInduced_Dipole(Eigen::Vector3d(-0.4, 0.3, 0.6));
    eeInteractor ee;
    const double e = ee.CalcPolarEnergy(seg1, seg2).E_indu_indu();
    const double e_ref = seg2[0].Induced_Charge() * seg2[0].Phi_noE() +
                         seg2[0].Induced_Dipole().dot(seg2[0].V_noE());
    BOOST_CHECK_CLOSE(e, e_ref, 1e-10);
    // and symmetric
    const double e21 = ee.CalcPolarEnergy(seg2, seg1).E_indu_indu();
    BOOST_CHECK_CLOSE(e, e21, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(thole_sites_unchanged) {
  // two Thole sites: FillDipoleInteraction is the Thole tensor
  PolarSegment a = OneSite(1, Eigen::Vector3d::Zero(), 9.0, false, 0.0);
  PolarSegment b = OneSite(2, Eigen::Vector3d(1.0, 2.0, 0.5), 6.0, false, 0.0);
  eeInteractor ee(0.39);
  BOOST_CHECK(ee.FillDipoleInteraction(a[0], b[0]) ==
              ee.FillTholeInteraction(a[0], b[0]));
  // Thole damping requested for charge-response dipoles
  PolarSegment c = OneSite(3, Eigen::Vector3d(1.0, 2.0, 0.5), 6.0, false, 0.0);
  c[0].setChargeResponse(1.0, true);
  BOOST_CHECK(ee.FillDipoleInteraction(a[0], c[0]) ==
              ee.FillTholeInteraction(a[0], c[0]));
  // no couplings inside a charge-response segment
  PolarSegment two("two", 4);
  two.push_back(PolarSite(0, "H", Eigen::Vector3d::Zero()));
  two.push_back(PolarSite(1, "H", Eigen::Vector3d(0, 0, 1.4)));
  for (PolarSite& s : two) {
    s.setpolarization(Eigen::Matrix3d::Identity());
    s.setInduced_Dipole(Eigen::Vector3d(0, 0, 1));
  }
  BOOST_CHECK(std::abs(ee.CalcPolarEnergy_IntraSegment(two)) > 1e-3);
  Eigen::MatrixXd pi(2, 2);
  pi << 1, -1, -1, 1;
  two.setChargeResponse(ChargeResponse({"H", "H"}, pi));
  BOOST_CHECK_EQUAL(ee.CalcPolarEnergy_IntraSegment(two), 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
