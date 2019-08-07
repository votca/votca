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
#include "votca/xtp/numerical_integrations.h"
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/polarregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/staticregion.h>

namespace votca {
namespace xtp {

void QMRegion::Initialize(const tools::Property& prop) {
  if (this->_id != 0) {
    throw std::runtime_error(
        this->identify() +
        " must always be region 0. Currently only one qm region is possible.");
  }

  std::string statestring =
      prop.ifExistsReturnElseThrowRuntimeError<std::string>("state");
  _initstate.FromString(statestring);
  if (_initstate.Type() == QMStateType::Hole ||
      _initstate.Type() == QMStateType::Electron) {
    throw std::runtime_error(
        "Charged QM Regions are not implemented currently");
  }
  if (_initstate.Type().isExciton() || _initstate.Type().isGWState()) {

    _do_gwbse = true;
    std::string gwbsexml =
        prop.ifExistsReturnElseThrowRuntimeError<std::string>("options_gwbse");
    tools::load_property_from_xml(_gwbseoptions, gwbsexml);
    if (prop.exists("filter")) {
      tools::Property filter = prop.get("filter");
      _filter.setLogger(&_log);
      _filter.Initialize(filter);
      _filter.setInitialState(_initstate);
      _filter.PrintInfo();
    } else {
      throw std::runtime_error("No filter for excited states specified");
    }
  }

  _grid_accuracy_for_ext_interaction = prop.ifExistsReturnElseReturnDefault(
      "grid_for_potential", _grid_accuracy_for_ext_interaction);
  _DeltaE = prop.ifExistsReturnElseReturnDefault("tolerance_energy", _DeltaE);
  _DeltaD = prop.ifExistsReturnElseReturnDefault("tolerance_density", _DeltaD);

  std::string dftxml =
      prop.ifExistsReturnElseThrowRuntimeError<std::string>("options_dft");
  tools::load_property_from_xml(_dftoptions, dftxml);
}

bool QMRegion::Converged() const {
  if (!_E_hist.filled()) {
    return false;
  }

  double Echange = _E_hist.getDiff();
  double Dchange =
      _Dmat_hist.getDiff().norm() / double(_Dmat_hist.back().cols());
  double Dmax = _Dmat_hist.getDiff().cwiseAbs().maxCoeff();
  std::string info = "not converged";
  bool converged = false;
  if (Dchange < _DeltaD && Dmax < _DeltaD && std::abs(Echange) < _DeltaE) {
    info = "converged";
    converged = true;
  }
  XTP_LOG_SAVE(logINFO, _log)
      << " Region:" << this->identify() << " " << this->getId() << " is "
      << info << " deltaE=" << Echange << " RMS Dmat=" << Dchange
      << " MaxDmat=" << Dmax << std::flush;
  return converged;
}

void QMRegion::Evaluate(std::vector<std::unique_ptr<Region> >& regions) {

  std::vector<double> interact_energies = ApplyInfluenceOfOtherRegions(regions);
  double e_ext =
      std::accumulate(interact_energies.begin(), interact_energies.end(), 0.0);
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp()
      << " Calculated interaction potentials with other regions. E[hrt]= "
      << e_ext << std::flush;
  XTP_LOG_SAVE(logINFO, _log) << "Writing inputs" << std::flush;
  _qmpackage->setRunDir(_workdir);
  _qmpackage->WriteInputFile(_orb);
  XTP_LOG_SAVE(logINFO, _log) << "Running DFT calculation" << std::flush;
  bool run_success = _qmpackage->Run();
  if (!run_success) {
    _info = false;
    _errormsg = "DFT-run failed";
    return;
  }

  bool Logfile_parse = _qmpackage->ParseLogFile(_orb);
  if (!Logfile_parse) {
    _info = false;
    _errormsg = "Parsing DFT logfile failed.";
    return;
  }
  bool Orbfile_parse = _qmpackage->ParseMOsFile(_orb);
  if (!Orbfile_parse) {
    _info = false;
    _errormsg = "Parsing DFT orbfile failed.";
    return;
  }
  QMState state = QMState("groundstate");
  double energy = _orb.getDFTTotalEnergy();
  if (_do_gwbse) {
    GWBSE gwbse(_orb);
    gwbse.setLogger(&_log);
    gwbse.Initialize(_gwbseoptions);
    gwbse.Evaluate();
    state = _filter.CalcStateAndUpdate(_orb);
    if (state.Type().isExciton()) {
      energy += _orb.getExcitedStateEnergy(state);
    } else {
      energy = _orb.getExcitedStateEnergy(state);
    }
  }
  _E_hist.push_back(energy);
  _Dmat_hist.push_back(_orb.DensityMatrixFull(state));
  return;
}

void QMRegion::push_back(const QMMolecule& mol) {
  if (_orb.QMAtoms().size() == 0) {
    _orb.QMAtoms() = mol;
  } else {
    _orb.QMAtoms().AddContainer(mol);
  }
  _size++;
}

double QMRegion::charge() const {
  double charge = 0.0;
  if (!_do_gwbse) {
    double nuccharge = 0.0;
    for (const QMAtom& a : _orb.QMAtoms()) {
      nuccharge += a.getNuccharge();
    }

    double electrons = _orb.getNumberOfAlphaElectrons() * 2;
    charge = nuccharge - electrons;
  } else {
    QMState state = _filter.InitialState();
    if (state.Type().isExciton()) {
      charge = 0.0;
    } else if (state.Type().isSingleParticleState()) {
      if (state.Index() <= _orb.getHomo()) {
        charge = +1.0;
      } else {
        charge = -1.0;
      }
    }
  }
  return charge;
}

void QMRegion::AppendResult(tools::Property& prop) const {
  prop.add("E_total", std::to_string(_E_hist.back() * tools::conv::hrt2ev));
  if (_do_gwbse) {
    prop.add("Initial_State", _filter.InitialState().ToString());
    prop.add("Final_State", _filter.CalcState(_orb).ToString());
  }
}

void QMRegion::Reset() {

  std::string dft_package_name =
      _dftoptions.get("package.name").as<std::string>();
  _qmpackage =
      std::unique_ptr<QMPackage>(QMPackages().Create(dft_package_name));
  _qmpackage->setLog(&_log);
  _qmpackage->Initialize(_dftoptions);
  int charge = 0;
  if (_initstate.Type() == QMStateType::Electron) {
    charge = -1;
  } else if (_initstate.Type() == QMStateType::Hole) {
    charge = +1;
  }
  _qmpackage->setCharge(charge);
  return;
}
double QMRegion::InteractwithQMRegion(const QMRegion& region) {
  throw std::runtime_error(
      "QMRegion-QMRegion interaction is not implemented yet.");
  return 0.0;
}
double QMRegion::InteractwithPolarRegion(const PolarRegion& region) {
  _qmpackage->AddRegion(region);
  return 0.0;
}
double QMRegion::InteractwithStaticRegion(const StaticRegion& region) {
  _qmpackage->AddRegion(region);
  return 0.0;
}

void QMRegion::WritePDB(csg::PDBWriter& writer) const {
  writer.WriteContainer(_orb.QMAtoms());
}

void QMRegion::AddNucleiFields(std::vector<PolarSegment>& segments,
                               const StaticSegment& seg) const {
  eeInteractor e;
#pragma omp parallel for
  for (int i = 0; i < int(segments.size()); ++i) {
    e.ApplyStaticField<StaticSegment, Estatic::noE_V>(seg, segments[i]);
  }
}

void QMRegion::ApplyQMFieldToPolarSegments(
    std::vector<PolarSegment>& segments) const {

  NumericalIntegration numint;
  AOBasis basis =
      _orb.SetupDftBasis();  // grid needs a basis in scope all the time
  numint.GridSetup(_grid_accuracy_for_ext_interaction, _orb.QMAtoms(), basis);

  QMState state = QMState("groundstate");
  if (_do_gwbse) {
    state = _filter.CalcState(_orb);
  }
  Eigen::MatrixXd dmat = _orb.DensityMatrixFull(state);
  double Ngrid = numint.IntegrateDensity(dmat);
  AOOverlap overlap;
  overlap.Fill(basis);
  double N_comp = dmat.cwiseProduct(overlap.Matrix()).sum();
  if (std::abs(Ngrid - N_comp) > 0.001) {
    XTP_LOG(logDEBUG, _log) << "=======================" << std::flush;
    XTP_LOG(logDEBUG, _log)
        << "WARNING: Calculated Densities at Numerical Grid, Number of "
           "electrons "
        << Ngrid << " is far away from the the real value " << N_comp
        << ", you should increase the accuracy of the integration grid."
        << std::flush;
  }
#pragma omp parallel for
  for (int i = 0; i < int(segments.size()); ++i) {
    PolarSegment& seg = segments[i];
    for (PolarSite& site : seg) {
      site.V_noE() += numint.IntegrateField(site.getPos());
    }
  }

  StaticSegment seg(_orb.QMAtoms().getName(), _orb.QMAtoms().getId());
  for (const QMAtom& atom : _orb.QMAtoms()) {
    seg.push_back(StaticSite(atom, atom.getNuccharge()));
  }
  AddNucleiFields(segments, seg);
}

void QMRegion::WriteToCpt(CheckpointWriter& w) const {
  w(_id, "id");
  w(identify(), "type");
  w(_size, "size");
  w(_do_gwbse, "GWBSE");
  w(_initstate.ToString(), "initial_state");
  w(_grid_accuracy_for_ext_interaction, "ext_grid");
  CheckpointWriter v = w.openChild("orbitals");
  _orb.WriteToCpt(v);

  CheckpointWriter v2 = w.openChild("E-hist");
  _E_hist.WriteToCpt(v2);

  CheckpointWriter v3 = w.openChild("D-hist");
  _Dmat_hist.WriteToCpt(v3);
  if (_do_gwbse) {
    CheckpointWriter v4 = w.openChild("statefilter");
    _filter.WriteToCpt(v4);
  }
}

void QMRegion::ReadFromCpt(CheckpointReader& r) {
  r(_id, "id");
  r(_size, "size");
  r(_do_gwbse, "GWBSE");
  std::string state;
  r(state, "initial_state");
  _initstate.FromString(state);
  r(_grid_accuracy_for_ext_interaction, "ext_grid");
  CheckpointReader rr = r.openChild("orbitals");
  _orb.ReadFromCpt(rr);

  CheckpointReader rr2 = r.openChild("E-hist");
  _E_hist.ReadFromCpt(rr2);

  CheckpointReader rr3 = r.openChild("D-hist");
  _Dmat_hist.ReadFromCpt(rr3);
  if (_do_gwbse) {
    CheckpointReader rr4 = r.openChild("statefilter");
    _filter.ReadFromCpt(rr4);
  }
}

}  // namespace xtp
}  // namespace votca
