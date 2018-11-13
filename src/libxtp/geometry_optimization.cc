/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/forces.h>
#include <votca/xtp/bfgs-trm.h>
#include <votca/xtp/energy_costfunction.h>
#include <functional>
#include <votca/xtp/statefilter.h>

namespace votca {
  namespace xtp {

    void GeometryOptimization::Initialize(tools::Property &options) {

      std::string statestring = options.ifExistsReturnElseThrowRuntimeError<std::string>(".state");
      _opt_state.FromString(statestring);
      if (!_opt_state.Type().isExciton()) {
        throw std::runtime_error("At the moment only excitonic states can be optimized");
      }
      // default convergence parameters from ORCA
      _conv.deltaE = options.ifExistsReturnElseReturnDefault<double>(".convergence.energy", 1.e-6); // Hartree
      _conv.RMSForce = options.ifExistsReturnElseReturnDefault<double>(".convergence.RMSForce", 3.e-5); // Hartree/Bohr
      _conv.MaxForce = options.ifExistsReturnElseReturnDefault<double>(".convergence.MaxForce", 1.e-4); // Hartree/Bohr
      _conv.RMSStep = options.ifExistsReturnElseReturnDefault<double>(".convergence.RMSStep", 6.e-4); // Bohr
      _conv.MaxStep = options.ifExistsReturnElseReturnDefault<double>(".convergence.MaxStep", 1.e-3); // Bohr
      _trust_radius = options.ifExistsReturnElseReturnDefault<double>(".trust", 0.01); // Angstrom
      _trust_radius *= tools::conv::ang2bohr; // initial trust radius in a.u.

      _max_iteration = options.ifExistsReturnElseReturnDefault<unsigned>(".maxiter", 50);

      _trajfile = options.ifExistsReturnElseReturnDefault<std::string>(".trajectory_file", "optimisation.trj");


      std::vector<std::string> choices = {"BFGS-TRM"};
      _optimizer = options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(".optimizer.method", choices);

      if (options.exists(".forces")) {
        _force_options = options.get(".forces");
      } else {
        throw std::runtime_error("No forces options provided");
      }
      if (options.exists(".filter")) {
        _filter_options = options.get(".filter");
      } else {
        throw std::runtime_error("No filter options set");
      }

      return;
    }

    void GeometryOptimization::Evaluate() {
      CTP_LOG(ctp::logINFO, *_pLog) << "Requested geometry optimization of excited state " << _opt_state.ToString() << std::flush;

      Statefilter filter;
      filter.Initialize(_filter_options);
      filter.setInitialState(_opt_state);
      filter.setLogger(_pLog);
      filter.PrintInfo();

      // get a force object
      Forces force_engine(_gwbse_engine, filter);
      force_engine.Initialize(_force_options);
      force_engine.setLog(_pLog);
      CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("Convergence of total energy: %1$8.6f Hartree ") % _conv.deltaE).str() << std::flush;
      CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("Convergence of RMS Force:    %1$8.6f Hartree/Bohr ") % _conv.RMSForce).str() << std::flush;
      CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("Convergence of Max Force:    %1$8.6f Hartree/Bohr ") % _conv.MaxForce).str() << std::flush;
      CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("Convergence of RMS Step:     %1$8.6f Bohr ") % _conv.RMSStep).str() << std::flush;
      CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("Convergence of Max Step:     %1$8.6f Bohr ") % _conv.MaxStep).str() << std::flush;
      CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("Initial trust radius:        %1$8.6f Bohr") % _trust_radius).str() << std::flush;

      Energy_costfunction e_cost = Energy_costfunction(_gwbse_engine, filter, _orbitals, force_engine);
      e_cost.setConvergenceParameters(_conv);
      e_cost.setLog(_pLog);
      // get the optimizer
      if (_optimizer == "BFGS-TRM") {
        BFGSTRM bfgstrm(e_cost);
        std::vector<std::function<void()> > callbacks;
        std::function<void()> reporting=std::bind(Report,std::cref(bfgstrm),std::cref(force_engine),std::cref(_pLog));
        callbacks.push_back(reporting);
        std::function<void()> filewrite=std::bind(WriteTrajectory,_trajfile, _orbitals.QMAtoms(), std::cref(bfgstrm));
        callbacks.push_back(filewrite);
        bfgstrm.setCallbacks(callbacks);
        bfgstrm.setNumofIterations(_max_iteration);
        bfgstrm.setTrustRadius(_trust_radius);
        bfgstrm.setLog(_pLog);
        bfgstrm.Optimize(Energy_costfunction::QMAtoms2Vector(_orbitals.QMAtoms()));
      }
      return;
    }

    void GeometryOptimization::Report(const BFGSTRM& bfgstrm,const Forces& forces,ctp::Logger* pLog){

      CTP_LOG(ctp::logINFO, *pLog) << std::flush;
      CTP_LOG(ctp::logINFO, *pLog) << (boost::format(" =========== OPTIMIZATION SUMMARY ================================= ")).str() << std::flush;
      CTP_LOG(ctp::logINFO, *pLog) << " At iteration  " << bfgstrm.getIteration() << std::flush;
      CTP_LOG(ctp::logINFO, *pLog) << (boost::format("   ---- POSITIONS (Angstrom)   ")).str() << std::flush;
      CTP_LOG(ctp::logINFO, *pLog) << (boost::format("   Atom\t x\t  y\t  z ")).str() << std::flush;
      const Eigen::VectorXd& atomvec = bfgstrm.getParameters();
      for (unsigned i = 0; i < atomvec.size(); i += 3) {
        CTP_LOG(ctp::logINFO, *pLog) << (boost::format(" %1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f")
                % (i/3) % (atomvec(i) * votca::tools::conv::bohr2ang) % (atomvec(i + 1) * votca::tools::conv::bohr2ang) % (atomvec(i + 2) * votca::tools::conv::bohr2ang)).str() << std::flush;
      }
      CTP_LOG(ctp::logINFO, *pLog) << (boost::format("   Total energy:     %1$12.8f Hartree ") % bfgstrm.getCost()).str() << std::flush;
      CTP_LOG(ctp::logINFO, *pLog) << (boost::format("   Trust radius:     %1$12.8f Bohr     ") % bfgstrm.getTrustRadius()).str() << std::flush;
      forces.Report();
      return;
    }

    

    void GeometryOptimization::WriteTrajectory(const std::string& filename, std::vector< QMAtom* >& atoms, const BFGSTRM& bfgstrm){
      std::ofstream ofs;
      if (bfgstrm.getIteration() == 0) {
        ofs.open(filename.c_str(), std::ofstream::out);
      } else {
        ofs.open(filename.c_str(), std::ofstream::app);
      }
      if (!ofs.is_open()) {
        throw std::runtime_error("Bad file handle: " + filename);
      }
      // write coordinates as xyz file
      ofs << atoms.size() << std::endl;
      ofs << "iteration " << bfgstrm.getIteration() << " energy " << bfgstrm.getCost() << " Hartree" << std::endl;
      Energy_costfunction::Vector2QMAtoms(bfgstrm.getParameters(), atoms);
      for (const QMAtom* atom : atoms) {
        tools::vec pos = atom->getPos() * tools::conv::bohr2ang;
        ofs << atom->getType() << " " << pos.getX() << " " << pos.getY() << " " << pos.getZ() << std::endl;
      }
      ofs.close();
      return;
    }

  }
}
