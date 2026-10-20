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

#define BOOST_TEST_MODULE chargeresponsesolver_test

// Standard includes
#include <fstream>
#include <memory>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/chargeresponsesolver.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/staticregion.h"

using namespace votca::xtp;
using namespace votca;

namespace {

// A reference case of the Python reference solver of the charge-response
// study (cr_reference.py, export_cpp): molecules, point charges, and the
// induced charges/dipoles and polarization energy they give.
struct RefCase {
  std::string dd;
  std::vector<PolarSegment> segments;
  double width = 0.0;  // bohr (the same for all CR molecules)
  std::vector<std::pair<Eigen::Vector3d, double>> charges;
  double E_pol = 0.0;
  Eigen::Matrix3d alpha = Eigen::Matrix3d::Zero();
  std::vector<std::pair<double, Eigen::Vector3d>> sites;
};

RefCase Load(const std::string& name) {
  std::ifstream in(std::string(XTP_TEST_DATA_FOLDER) + "/chargeresponse/" +
                   name + ".dat");
  BOOST_REQUIRE(in.is_open());
  RefCase c;
  std::string key;
  Index n = 0;
  in >> key >> c.dd;
  in >> key >> n;
  for (Index m = 0; m < n; ++m) {
    std::string kind;
    Index na = 0;
    double width = 0.0;
    in >> kind >> na >> width;
    PolarSegment seg("mol" + std::to_string(m), m);
    std::vector<std::string> elements;
    for (Index i = 0; i < na; ++i) {
      std::string el;
      Eigen::Vector3d pos;
      double xx, xy, xz, yy, yz, zz;
      in >> el >> pos(0) >> pos(1) >> pos(2) >> xx >> xy >> xz >> yy >> yz >>
          zz;
      PolarSite site(i, el, pos);
      Eigen::Matrix3d a;
      a << xx, xy, xz, xy, yy, yz, xz, yz, zz;
      site.setpolarization(a);
      seg.push_back(site);
      elements.push_back(el);
    }
    if (kind == "cr") {
      Eigen::MatrixXd pi(na, na);
      for (Index i = 0; i < na; ++i) {
        for (Index j = 0; j < na; ++j) {
          in >> pi(i, j);
        }
      }
      seg.setChargeResponse(ChargeResponse(elements, pi));
      c.width = width;
    }
    c.segments.push_back(seg);
  }
  in >> key >> n;
  for (Index k = 0; k < n; ++k) {
    Eigen::Vector3d p;
    double q;
    in >> p(0) >> p(1) >> p(2) >> q;
    c.charges.push_back({p, q});
  }
  in >> key >> c.E_pol;
  in >> key;
  for (Index k = 0; k < 9; ++k) {
    in >> c.alpha(k / 3, k % 3);
  }
  in >> key >> n;
  for (Index k = 0; k < n; ++k) {
    double dq;
    Eigen::Vector3d mu;
    in >> dq >> mu(0) >> mu(1) >> mu(2);
    c.sites.push_back({dq, mu});
  }
  BOOST_REQUIRE(in.good());
  return c;
}

tools::Property PolarOptions(const RefCase& c) {
  tools::Property prop;
  prop.add("max_iter", "500");
  prop.add("tolerance_dipole", "1e-12");
  prop.add("tolerance_energy", "1e-10");
  prop.add("exp_damp", "0.39");
  std::ostringstream w;
  w.precision(17);
  w << c.width / tools::conv::ang2bohr;
  prop.add("cr_width", w.str());
  prop.add("cr_dipole_damping", c.dd);
  return prop;
}

// Runs the case as a polar region (the molecules) next to a static region
// (the point charges) and compares with the reference.
void CheckCase(const std::string& name) {
  RefCase c = Load(name);
  Logger log;
  auto polar = std::make_unique<PolarRegion>(0, log);
  for (const PolarSegment& seg : c.segments) {
    polar->push_back(seg);
  }
  polar->Initialize(PolarOptions(c));
  auto stat = std::make_unique<StaticRegion>(1, log);
  StaticSegment charges("charges", 100);
  for (Index k = 0; k < Index(c.charges.size()); ++k) {
    StaticSite site(k, "H", c.charges[k].first);
    site.setCharge(c.charges[k].second);
    charges.push_back(site);
  }
  stat->push_back(charges);
  std::vector<std::unique_ptr<Region>> regions;
  regions.push_back(std::move(polar));
  regions.push_back(std::move(stat));
  regions[0]->Reset();
  regions[0]->Evaluate(regions);
  const PolarRegion& p = dynamic_cast<const PolarRegion&>(*regions[0]);
  BOOST_REQUIRE(p.Successful());
  // the molecules have no permanent multipoles: the region's energy is the
  // polarization energy 1/2 x^T b
  BOOST_CHECK_CLOSE(p.Etotal(), c.E_pol, 1e-6);
  Index k = 0;
  double max_dq = 0.0, max_dmu = 0.0, ref_q = 0.0, ref_mu = 0.0;
  for (const PolarSegment& seg : p) {
    for (const PolarSite& site : seg) {
      max_dq =
          std::max(max_dq, std::abs(site.Induced_Charge() - c.sites[k].first));
      max_dmu =
          std::max(max_dmu, (site.Induced_Dipole() - c.sites[k].second).norm());
      ref_q = std::max(ref_q, std::abs(c.sites[k].first));
      ref_mu = std::max(ref_mu, c.sites[k].second.norm());
      ++k;
    }
  }
  BOOST_CHECK_EQUAL(k, Index(c.sites.size()));
  BOOST_CHECK_SMALL(max_dq, 1e-7 * std::max(ref_q, 1e-3));
  BOOST_CHECK_SMALL(max_dmu, 1e-7 * ref_mu);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(chargeresponsesolver_test)

BOOST_AUTO_TEST_CASE(single_molecule_polarizability) {
  // a uniform field on one charge-response molecule: the induced dipole is
  // the fitted polarizability X^T Pi X + sum alpha_res
  RefCase c = Load("pentacene_single_gaussian");
  std::vector<PolarSegment> segs = c.segments;
  for (PolarSite& site : segs[0]) {
    site.setChargeResponse(c.width, false);
  }
  eeInteractor ee(0.39);
  Eigen::Matrix3d alpha;
  for (Index d = 0; d < 3; ++d) {
    const Eigen::Vector3d F = Eigen::Vector3d::Unit(d);
    for (PolarSite& site : segs[0]) {
      site.Reset();
      site.V() = -F;                       // grad phi
      site.Phi() = -F.dot(site.getPos());  // phi = -F.r
    }
    ChargeResponseSolver solver(segs, ee);
    BOOST_REQUIRE(solver.Solve(1e-12, 100, false).converged);
    Eigen::Vector3d mu = Eigen::Vector3d::Zero();
    for (const PolarSite& site : segs[0]) {
      mu += site.Induced_Charge() * site.getPos() + site.Induced_Dipole();
    }
    alpha.col(d) = mu;
  }
  BOOST_CHECK(alpha.isApprox(c.alpha, 1e-8));
}

BOOST_AUTO_TEST_CASE(operator_symmetric_positive) {
  RefCase c = Load("water_trimer_gaussian");
  std::vector<PolarSegment> segs = c.segments;
  for (PolarSegment& seg : segs) {
    for (PolarSite& site : seg) {
      site.setChargeResponse(c.width, false);
    }
  }
  eeInteractor ee(0.39);
  ChargeResponseSolver solver(segs, ee);
  const Index n = solver.size();
  Eigen::MatrixXd A(n, n);
  for (Index i = 0; i < n; ++i) {
    A.col(i) = solver.Apply(Eigen::VectorXd::Unit(n, i));
  }
  BOOST_CHECK_SMALL((A - A.transpose()).cwiseAbs().maxCoeff(), 1e-12);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  BOOST_CHECK(es.eigenvalues()(0) > 0.5);
}

BOOST_AUTO_TEST_CASE(reference_cases) {
  CheckCase("pentacene_dimer_pair0_gaussian");
  CheckCase("pentacene_mixed_CR_A_Thole_B_thole");
  CheckCase("water_trimer_gaussian");
  CheckCase("water_trimer_thole");
}

BOOST_AUTO_TEST_CASE(induced_charges_checkpoint) {
  RefCase c = Load("water_trimer_gaussian");
  PolarSegment seg = c.segments[0];
  for (Index i = 0; i < seg.size(); ++i) {
    seg[i].setInduced_Charge(0.1 * double(i) - 0.1);
  }
  {
    CheckpointFile ff("crsolver.hdf5");
    CheckpointWriter ww = ff.getWriter();
    seg.WriteToCpt(ww);
  }
  CheckpointFile ff("crsolver.hdf5", CheckpointAccessLevel::READ);
  CheckpointReader rr = ff.getReader();
  PolarSegment seg2(rr);
  for (Index i = 0; i < seg.size(); ++i) {
    BOOST_CHECK_EQUAL(seg2[i].Induced_Charge(), seg[i].Induced_Charge());
  }
}

BOOST_AUTO_TEST_SUITE_END()
