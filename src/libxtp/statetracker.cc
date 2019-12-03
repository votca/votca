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

#include "votca/xtp/aomatrix.h"
#include <numeric>
#include <votca/xtp/statetracker.h>

namespace votca {
namespace xtp {
using std::flush;

void StateTracker::Initialize(tools::Property& options) {
  if (options.exists("oscillator_strength")) {
    _use_osctracker = true;
    _oscthreshold = options.ifExistsReturnElseThrowRuntimeError<double>(
        "oscillator_strength");
  }
  if (options.exists("overlap")) {
    _use_overlaptracker = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("overlap", 0.0);
  }
  if (options.exists("localization")) {
    _use_localizationtracker = true;
    std::string indices =
        options.ifExistsReturnElseThrowRuntimeError<std::string>(
            "localization.fragment");
    _fragment_loc = QMFragment<double>(0, indices);
    _fragment_loc.value() = options.ifExistsReturnElseThrowRuntimeError<double>(
        "localization.threshold");
  }

  if (options.exists("charge_transfer")) {
    _use_dQtracker = true;
    std::string indices =
        options.ifExistsReturnElseThrowRuntimeError<std::string>(
            "charge_transfer.fragment");
    _fragment_dQ = QMFragment<double>(0, indices);
    _fragment_dQ.value() = options.ifExistsReturnElseThrowRuntimeError<double>(
        "charge_transfer.threshold");
  }
  if (_use_dQtracker && _use_localizationtracker) {
    throw std::runtime_error(
        "Cannot use localization and charge_transfer tracker at the same "
        "time.");
  }
  if (options.exists("overlapbse")) {
    _use_overlaptracker_bse = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("overlapbse", 0.0);
  }
  if (options.exists("density")) {
    _use_densitytracker = true;
    _dmatthreshold =
        options.ifExistsReturnElseReturnDefault<double>("density", 0.0);
  }
}

void StateTracker::PrintInfo() const {
  XTP_LOG(Log::error, *_log)
      << "Initial state: " << _statehist[0].ToString() << flush;
  if (_statehist.size() > 1) {
    XTP_LOG(Log::error, *_log)
        << "Last state: " << _statehist.back().ToString() << flush;
  }
  if (_use_osctracker) {
    XTP_LOG(Log::error, *_log)
        << "Using oscillator strength tracker with threshold " << _oscthreshold
        << flush;
  }
  if (_use_overlaptracker) {
    if (_overlapthreshold == 0.0) {
      XTP_LOG(Log::error, *_log)
          << "Using overlap filer with no threshold " << flush;
    } else {
      XTP_LOG(Log::error, *_log) << "Using overlap filer with threshold "
                                 << _overlapthreshold << flush;
    }
  }
  if (_use_localizationtracker) {
    XTP_LOG(Log::error, *_log)
        << "Using localization tracker for " << _fragment_loc << flush;
  }
  if (_use_dQtracker) {
    XTP_LOG(Log::error, *_log)
        << "Using Delta Q tracker for fragment " << _fragment_dQ << flush;
  }
  if (_use_osctracker && _use_dQtracker) {
    XTP_LOG(Log::error, *_log) << "WARNING: trackering for optically active CT "
                                  "transition - might not make sense... "
                               << flush;
  }
  if (_use_overlaptracker_bse) {
    if (_overlapthreshold == 0.0) {
      XTP_LOG(Log::error, *_log)
          << "Using overlap tracker BSE with no threshold " << flush;
    } else {
      XTP_LOG(Log::error, *_log) << "Using overlap tracker BSE with threshold "
                                 << _overlapthreshold << flush;
    }
  }
  if (_use_densitytracker) {
    if (_dmatthreshold == 0.0) {
      XTP_LOG(Log::error, *_log)
          << "Using density filter with no threshold " << flush;
    } else {
      XTP_LOG(Log::error, *_log)
          << "Using density filter with threshold " << _dmatthreshold << flush;
    }
  }

  if (_use_dQtracker + _use_osctracker + _use_localizationtracker +
          _use_overlaptracker + _use_overlaptracker_bse + _use_densitytracker <
      1) {
    XTP_LOG(Log::error, *_log) << "WARNING: No tracker is used " << flush;
  }
}

std::vector<Index> StateTracker::ComparePairofVectors(
    std::vector<Index>& vec1, std::vector<Index>& vec2) const {
  std::vector<Index> result(std::min(vec1, vec2));
  std::vector<Index>::iterator it;
  std::sort(vec1.begin(), vec1.end());
  std::sort(vec2.begin(), vec2.end());
  it = std::set_intersection(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                             result.begin());
  result.resize(it - result.begin());
  return result;
}


std::vector<Index> StateTracker::CollapseResults(
    std::vector<std::vector<Index> >& results) const {
  if (results.size() == 1) {
    return results[0];
  } else {
    std::vector<Index> result = results[0];
    for (Index i = 1; i < Index(results.size()); i++) {
      result = ComparePairofVectors(result, results[i]);
    }
    return result;
  }
}

QMState StateTracker::CalcState(const Orbitals& orbitals) const {

  if (_use_dQtracker + _use_overlaptracker + _use_localizationtracker +
          _use_osctracker <
      1) {
    return _statehist[0];
  }

  std::vector<std::vector<Index> > results;
  if (_use_osctracker) {
    results.push_back(OscTracker(orbitals));
  }
  if (_use_localizationtracker) {
    results.push_back(LocTracker(orbitals));
  }
  if (_use_overlaptracker) {
    results.push_back(OverlapTracker(orbitals));
  }
  if (_use_dQtracker) {
    results.push_back(DeltaQTracker(orbitals));
  }
  if (_use_overlaptracker_bse) {
    results.push_back(OverlapTrackerBSE(orbitals));
  }
  if (_use_densitytracker) {
    results.push_back(DensityTracker(orbitals));
  }

  std::vector<Index> result = CollapseResults(results);
  QMState state;
  if (result.size() < 1) {
    state = _statehist.back();
    XTP_LOG(Log::error, *_log)
        << "No State found by tracker using last state: " << state.ToString()
        << flush;
  } else {
    state = QMState(_statehist.back().Type(), result[0], false);
    XTP_LOG(Log::error, *_log)
        << "Next State is: " << state.ToString() << flush;
  }
  return state;
}

QMState StateTracker::CalcStateAndUpdate(const Orbitals& orbitals) {
  QMState result = CalcState(orbitals);
  _statehist.push_back(result);
  if (_use_overlaptracker) {
    UpdateLastCoeff(orbitals);
  }
  if (_use_overlaptracker_bse) {
    UpdateLastBSE_R(orbitals);
    UpdateLastCoeff_matrix(orbitals);
  }
  if (_use_densitytracker) {
    UpdateLastDmat(orbitals);
  }
  return result;
}

std::vector<Index> StateTracker::OscTracker(const Orbitals& orbitals) const {
  Eigen::VectorXd oscs = orbitals.Oscillatorstrengths();
  std::vector<Index> indexes;
  for (Index i = 0; i < oscs.size(); i++) {
    if (oscs[i] > _oscthreshold) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

std::vector<Index> StateTracker::LocTracker(const Orbitals& orbitals) const {
  std::vector<Index> indexes;
  Lowdin low;
  QMFragment<BSE_Population> frag;
  frag.copy_withoutvalue(_fragment_loc);
  std::vector<QMFragment<BSE_Population> > loc = {frag};
  low.CalcChargeperFragment(loc, orbitals, _statehist[0].Type());
  const Eigen::VectorXd& popE = loc[0].value().E;
  const Eigen::VectorXd& popH = loc[0].value().H;
  for (Index i = 0; i < popE.size(); i++) {
    if (std::abs(popE[i]) > _fragment_loc.value() &&
        std::abs(popH[i]) > _fragment_loc.value()) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

std::vector<Index> StateTracker::DeltaQTracker(const Orbitals& orbitals) const {
  std::vector<Index> indexes;
  Lowdin low;
  QMFragment<BSE_Population> frag;
  frag.copy_withoutvalue(_fragment_dQ);
  std::vector<QMFragment<BSE_Population> > loc = {frag};
  low.CalcChargeperFragment(loc, orbitals, _statehist[0].Type());
  Eigen::VectorXd dq = (loc[0].value().H + loc[0].value().E).cwiseAbs();

  for (Index i = 0; i < dq.size(); i++) {
    if (dq[i] > _fragment_dQ.value()) {
      indexes.push_back(i);
    }
  }
  return indexes;
}
Eigen::VectorXd StateTracker::CalculateOverlapBSE(
    const Orbitals& orbitals) const {
  // Define some useful Indeegers as the numeber of states, the number of occ and
  // virt level
  Index nostates = orbitals.NumberofStates(_statehist[0].Type());
  Index v = orbitals.getBSEvmax() + 1;
  Index c = orbitals.getBSEcmax() + 1;
  // Degine the AO overlap which is not step dependent unless you do geo.
  // optimization
  AOOverlap S_ao;
  S_ao.Fill(orbitals.SetupDftBasis());
  // Get the whole MO matrix at this step (n)
  Eigen::MatrixXd MO_mat = orbitals.MOs().eigenvectors();
  // Compute C as MO_mat(step n-1) * S_ao * MO_mat(step n)
  Eigen::MatrixXd C = _laststatecoeff_mat.transpose() * S_ao.Matrix() * MO_mat;

  // Get the diagonal blocks (to be updated with off diagonal block and beyond
  // TDA overlap
  Eigen::MatrixXd M_vv = C.block(0, 0, v, v);
  Eigen::MatrixXd M_cc = C.block(v, v, c - v, c - v);
  // Let's define the Reference vector containing the singlet amplitudes at step
  // n-1
  Eigen::VectorXd Ref = _lastbse_R;
  // I have to reshape this vector in a matrix form
  Eigen::MatrixXd A_vc_old = Eigen::Map<Eigen::MatrixXd>(Ref.data(), v, c - v);
  // Compute A_vc * M_cc'
  Eigen::MatrixXd two = A_vc_old * M_cc;

  // Inizialise the BSE overlap vector
  Eigen::VectorXd overlap_bse = Eigen::VectorXd::Zero(nostates);
  for (Index i = 0; i < nostates; i++) {
    QMState state(_statehist[0].Type(), i, false);
    Eigen::VectorXd a =
        orbitals.BSESinglets().eigenvectors().col(state.StateIdx());
    Eigen::MatrixXd A_vc = Eigen::Map<Eigen::MatrixXd>(a.data(), v, c - v);

    Eigen::MatrixXd one = M_vv * A_vc;

    Eigen::MatrixXd O = two.transpose() * one;
    double ov = (O).trace();
    overlap_bse(i) = ov;
  }
  return overlap_bse.cwiseAbs();
}

Eigen::VectorXd StateTracker::CalculateOverlap(const Orbitals& orbitals) const {
  Eigen::MatrixXd coeffs = CalcOrthoCoeffs(orbitals);
  Eigen::VectorXd overlap = (coeffs * _laststatecoeff).cwiseAbs2();
  return overlap;
}

Eigen::VectorXd StateTracker::CalculateDNorm(const Orbitals& orbitals) const {
  Index nostates = orbitals.NumberofStates(_statehist[0].Type());
  Eigen::VectorXd norm = Eigen::VectorXd::Zero(nostates);
  for (Index i = 0; i < nostates; i++) {
    QMState state(_statehist[0].Type(), i, false);
    Eigen::MatrixXd diff = (orbitals.DensityMatrixFull(state) - _lastdmat);
    norm(i) = diff.norm() / (_lastdmat.norm());
  }
  return norm;
}

Eigen::MatrixXd StateTracker::CalcOrthoCoeffs(const Orbitals& orbitals) const {
  QMStateType type = _statehist[0].Type();
  Eigen::MatrixXd coeffs;
  if (type.isSingleParticleState()) {
    if (type == QMStateType::DQPstate) {
      coeffs = orbitals.CalculateQParticleAORepresentation();
    } else {
      coeffs = orbitals.MOs().eigenvectors();
    }
  } else {
    throw std::runtime_error(
        "Overlap for excitons is implemented within TDA. Use a different "
        "filter");
  }
  return coeffs;
}

void StateTracker::UpdateLastCoeff(const Orbitals& orbitals) {
  Eigen::MatrixXd ortho_coeffs = CalcOrthoCoeffs(orbitals);
  Index offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }
  _laststatecoeff = ortho_coeffs.col(_statehist.back().StateIdx() - offset);
}
void StateTracker::UpdateLastCoeff_matrix(const Orbitals& orbitals) {
  _laststatecoeff_mat = orbitals.MOs().eigenvectors();
}

void StateTracker::UpdateLastDmat(const Orbitals& orbitals) {
  _lastdmat = orbitals.DensityMatrixFull(_statehist.back());
}

void StateTracker::UpdateLastBSE_R(const Orbitals& orbitals) {
  _lastbse_R =
      orbitals.BSESinglets().eigenvectors().col(_statehist.back().StateIdx());
}

std::vector<Index> StateTracker::OverlapTracker(
    const Orbitals& orbitals) const {
  std::vector<Index> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<Index>{_statehist[0].StateIdx()};
    return indexes;
  }

  Eigen::VectorXd Overlap = CalculateOverlap(orbitals);
  Index validelements = Index(Overlap.size());
  for (Index i = 0; i < Index(Overlap.size()); i++) {
    if (Overlap(i) < _overlapthreshold) {
      validelements--;
    }
  }

  std::vector<Index> index = std::vector<Index>(Overlap.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&Overlap](Index i1, Index i2) {
    return Overlap[i1] > Overlap[i2];
  });

  Index offset = 0;
  if (_statehist[0].Type().isGWState()) {
    offset = orbitals.getGWAmin();
  }

  for (Index i : index) {
    if (Index(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}
std::vector<Index> StateTracker::OverlapTrackerBSE(
    const Orbitals& orbitals) const {
  std::vector<Index> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<Index>{_statehist[0].StateIdx()};
    return indexes;
  }

  Eigen::VectorXd Overlap = CalculateOverlapBSE(orbitals);
  Index validelements = Overlap.size();
  for (Index i = 0; i < Overlap.size(); i++) {
    if (Overlap(i) < _overlapthreshold) {
      validelements--;
    }
  }

  std::vector<Index> index = std::vector<Index>(Overlap.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&Overlap](Index i1, Index i2) {
    return Overlap[i1] > Overlap[i2];
  });

  Index offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (Index i : index) {
    if (Index(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}

std::vector<Index> StateTracker::DensityTracker(
    const Orbitals& orbitals) const {
  std::vector<Index> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<Index>{_statehist[0].StateIdx()};
    return indexes;
  }

  Eigen::VectorXd dnorm = CalculateDNorm(orbitals);

  Index validelements = dnorm.size();
  for (Index i = 0; i < dnorm.size(); i++) {
    if (dnorm(i) > _dmatthreshold) {
      validelements--;
    }
  }
  std::vector<Index> index = std::vector<Index>(dnorm.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(),
                   [&dnorm](Index i1, Index i2) { return dnorm[i1] < dnorm[i2]; });
  Index offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (Index i : index) {
    if (Index(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}
void StateTracker::WriteToCpt(CheckpointWriter& w) const {
  std::vector<std::string> statehiststring;
  statehiststring.reserve(_statehist.size());
  for (const QMState& s : _statehist) {
    statehiststring.push_back(s.ToString());
  }
  w(statehiststring, "statehist");
  w(_use_osctracker, "osctracker");
  w(_oscthreshold, "oscthreshold");

  w(_use_overlaptracker_bse, "overlaptrackerbse");
  w(_lastbse_R, "lastbse_R");
  w(_use_densitytracker, "densitytracker");
  w(_dmatthreshold, "dmatthreshold");
  w(_lastdmat, "lastdmat");

  w(_use_overlaptracker, "overlaptracker");
  w(_overlapthreshold, "overlapthreshold");
  w(_laststatecoeff, "laststatecoeff");

  w(_use_localizationtracker, "localizationtracker");
  CheckpointWriter ww = w.openChild("fragment_loc");
  _fragment_loc.WriteToCpt(ww);

  w(_use_dQtracker, "dQtracker");
  CheckpointWriter ww2 = w.openChild("fragment_dQ");
  _fragment_dQ.WriteToCpt(ww2);
}

void StateTracker::ReadFromCpt(CheckpointReader& r) {
  std::vector<std::string> statehiststring;
  r(statehiststring, "statehist");
  _statehist.clear();
  _statehist.reserve(statehiststring.size());
  for (const std::string& s : statehiststring) {
    _statehist.push_back(QMState(s));
  }
  r(_use_osctracker, "osctracker");
  r(_oscthreshold, "oscthreshold");

  r(_use_overlaptracker_bse, "overlaptrackerbse");
  r(_lastbse_R, "lastbse_R");
  r(_use_densitytracker, "densitytracker");
  r(_dmatthreshold, "dmatthreshold");
  r(_lastdmat, "lastdmat");

  r(_use_overlaptracker, "overlaptracker");
  r(_overlapthreshold, "overlapthreshold");
  r(_laststatecoeff, "laststatecoeff");

  r(_use_localizationtracker, "localizationtracker");

  CheckpointReader rr = r.openChild("fragment_loc");
  _fragment_loc.ReadFromCpt(rr);

  r(_use_dQtracker, "dQtracker");
  CheckpointReader rr2 = r.openChild("fragment_dQ");
  _fragment_dQ.ReadFromCpt(rr2);
}

}  // namespace xtp
}  // namespace votca
