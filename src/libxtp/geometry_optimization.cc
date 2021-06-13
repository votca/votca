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
#include "functional"

// Local VOTCA includes
#include "votca/xtp/bfgs_trm.h"
#include "votca/xtp/energy_costfunction.h"
#include "votca/xtp/forces.h"
#include "votca/xtp/geometry_optimization.h"
#include "votca/xtp/statetracker.h"

namespace votca {
namespace xtp {

void GeometryOptimization::Initialize(tools::Property& options) {

  opt_state_ = options.get(".state").as<QMState>();
  if (!opt_state_.Type().isExciton()) {
    throw std::runtime_error(
        "At the moment only excitonic states can be optimized");
  }
  // default convergence parameters from ORCA
  conv_.deltaE = options.get(".convergence.energy").as<double>();  // Hartree
  conv_.RMSForce =
      options.get(".convergence.RMSForce").as<double>();  // Hartree/Bohr
  conv_.MaxForce =
      options.get(".convergence.MaxForce").as<double>();  // Hartree/Bohr
  conv_.RMSStep = options.get(".convergence.RMSStep").as<double>();  // Bohr
  conv_.MaxStep = options.get(".convergence.MaxStep").as<double>();  // Bohr
  trust_radius_ = options.get("optimizer.trust").as<double>();       // Angstrom
  trust_radius_ *= tools::conv::ang2bohr;  // initial trust radius in a.u.

  max_iteration_ = options.get(".maxiter").as<Index>();
  trajfile_ = options.get(".trajectory_file").as<std::string>();
  optimizer_ = options.get(".optimizer.method").as<std::string>();
  force_options_ = options.get(".forces");
  statetracker_options_ = options.get(".statetracker");
}

void GeometryOptimization::Evaluate() {
  XTP_LOG(Log::error, *pLog_)
      << "Requested geometry optimization of excited state "
      << opt_state_.ToString() << std::flush;

  StateTracker tracker;
  tracker.Initialize(statetracker_options_);
  tracker.setInitialState(opt_state_);
  tracker.setLogger(pLog_);
  tracker.PrintInfo();

  // get a force object
  Forces force_engine(gwbse_engine_, tracker);
  force_engine.Initialize(force_options_);
  force_engine.setLog(pLog_);
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("Convergence of total energy: %1$8.6f Hartree ") %
          conv_.deltaE)
             .str()
      << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("Convergence of RMS Force:    %1$8.6f Hartree/Bohr ") %
          conv_.RMSForce)
             .str()
      << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("Convergence of Max Force:    %1$8.6f Hartree/Bohr ") %
          conv_.MaxForce)
             .str()
      << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("Convergence of RMS Step:     %1$8.6f Bohr ") %
          conv_.RMSStep)
             .str()
      << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("Convergence of Max Step:     %1$8.6f Bohr ") %
          conv_.MaxStep)
             .str()
      << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("Initial trust radius:        %1$8.6f Bohr") %
          trust_radius_)
             .str()
      << std::flush;

  Energy_costfunction e_cost =
      Energy_costfunction(gwbse_engine_, tracker, orbitals_, force_engine);
  e_cost.setConvergenceParameters(conv_);
  e_cost.setLog(pLog_);
  // get the optimizer
  if (optimizer_ == "BFGS-TRM") {
    BFGSTRM bfgstrm(e_cost);
    std::vector<std::function<void()> > callbacks;
    std::function<void()> reporting = std::bind(
        Report, std::cref(bfgstrm), std::cref(force_engine), std::ref(*pLog_));
    callbacks.push_back(reporting);
    std::function<void()> filewrite = std::bind(
        WriteTrajectory, trajfile_, orbitals_.QMAtoms(), std::cref(bfgstrm));
    callbacks.push_back(filewrite);
    bfgstrm.setCallbacks(callbacks);
    bfgstrm.setNumofIterations(max_iteration_);
    bfgstrm.setTrustRadius(trust_radius_);
    bfgstrm.setLog(pLog_);
    bfgstrm.Optimize(Energy_costfunction::QMAtoms2Vector(orbitals_.QMAtoms()));
  }
  return;
}

void GeometryOptimization::Report(const BFGSTRM& bfgstrm, const Forces& forces,
                                  Logger& pLog) {

  XTP_LOG(Log::error, pLog) << std::flush;
  XTP_LOG(Log::error, pLog)
      << (boost::format("=========== OPTIMIZATION SUMMARY "
                        "================================= "))
             .str()
      << std::flush;
  XTP_LOG(Log::error, pLog)
      << "At iteration  " << bfgstrm.getIteration() << std::flush;
  XTP_LOG(Log::error, pLog)
      << (boost::format(" ---- POSITIONS (Angstrom)   ")).str() << std::flush;
  XTP_LOG(Log::error, pLog)
      << (boost::format(" Atom\t x\t  y\t  z ")).str() << std::flush;
  const Eigen::VectorXd& atomvec = bfgstrm.getParameters();
  for (Index i = 0; i < atomvec.size(); i += 3) {
    XTP_LOG(Log::error, pLog)
        << (boost::format("%1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f") % (i / 3) %
            (atomvec(i) * votca::tools::conv::bohr2ang) %
            (atomvec(i + 1) * votca::tools::conv::bohr2ang) %
            (atomvec(i + 2) * votca::tools::conv::bohr2ang))
               .str()
        << std::flush;
  }
  XTP_LOG(Log::error, pLog)
      << (boost::format("   Total energy:     %1$12.8f Hartree ") %
          bfgstrm.getCost())
             .str()
      << std::flush;
  XTP_LOG(Log::error, pLog)
      << (boost::format("   Trust radius:     %1$12.8f Bohr     ") %
          bfgstrm.getTrustRadius())
             .str()
      << std::flush;
  forces.Report();
  return;
}

void GeometryOptimization::WriteTrajectory(const std::string& filename,
                                           QMMolecule& atoms,
                                           const BFGSTRM& bfgstrm) {
  std::ofstream ofs;
  if (bfgstrm.getIteration() == 0) {
    ofs.open(filename, std::ofstream::out);
  } else {
    ofs.open(filename, std::ofstream::app);
  }
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + filename);
  }
  // write coordinates as xyz file
  ofs << atoms.size() << std::endl;
  ofs << "iteration " << bfgstrm.getIteration() << " energy "
      << bfgstrm.getCost() << " Hartree" << std::endl;
  Energy_costfunction::Vector2QMAtoms(bfgstrm.getParameters(), atoms);
  for (const QMAtom& atom : atoms) {
    const Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    ofs << atom.getElement() << " " << pos.x() << " " << pos.y() << " "
        << pos.z() << std::endl;
  }
  ofs.close();
  return;
}

}  // namespace xtp
}  // namespace votca
