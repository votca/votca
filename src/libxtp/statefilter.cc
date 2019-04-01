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
#include <votca/xtp/statefilter.h>

namespace votca {
namespace xtp {
using std::flush;

void Statefilter::Initialize(tools::Property& options) {
  if (options.exists("oscillator_strength")) {
    _use_oscfilter = true;
    _oscthreshold = options.ifExistsReturnElseThrowRuntimeError<double>(
        "oscillator_strength");
  }
  if (options.exists("overlap")) {
    _use_overlapfilter = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("overlap", 0.0);
  }
  if (options.exists("localisation")) {
    _use_localisationfilter = true;
    std::string temp = options.get("localisation").as<std::string>();
    tools::Tokenizer tok_cleanup(temp, ", \n\t");
    std::vector<std::string> strings_vec;
    tok_cleanup.ToVector(strings_vec);
    if (strings_vec.size() != 2) {
      throw std::runtime_error(
          "statefiler: Fragment and localisation threshold are not separated");
    }
    if (strings_vec[0] == "a" || strings_vec[0] == "A") {
      _localiseonA = true;
    } else if (strings_vec[0] == "b" || strings_vec[0] == "B") {
      _localiseonA = false;
    } else {
      throw std::runtime_error(
          "statefiler: Fragment label not known, either A or B");
    }
    _loc_threshold = boost::lexical_cast<double>(strings_vec[1]);
  }

  if (options.exists("charge_transfer")) {
    _use_dQfilter = true;
    _dQ_threshold =
        options.ifExistsReturnElseThrowRuntimeError<double>("charge_transfer");
  }
  if (_use_dQfilter && _use_localisationfilter) {
    throw std::runtime_error(
        "Cannot use localisation and charge_transfer filter at the same time.");
  }
}

void Statefilter::PrintInfo() const {
  CTP_LOG(ctp::logDEBUG, *_log)
      << "Initial state: " << _statehist[0].ToString() << flush;
  if (_statehist.size() > 1) {
    CTP_LOG(ctp::logDEBUG, *_log)
        << "Last state: " << _statehist.back().ToString() << flush;
  }
  if (_use_oscfilter) {
    CTP_LOG(ctp::logDEBUG, *_log)
        << "Using oscillator strength filter with cutoff " << _oscthreshold
        << flush;
  }
  if (_use_overlapfilter) {
    if (_overlapthreshold == 0.0) {
      CTP_LOG(ctp::logDEBUG, *_log)
          << "Using overlap filer with no cutoff " << flush;
    } else {
      CTP_LOG(ctp::logDEBUG, *_log)
          << "Using overlap filer with cutoff " << _overlapthreshold << flush;
    }
  }
  if (_use_localisationfilter) {
    std::string fragment = "A";
    if (!_localiseonA) {
      fragment = "B";
    }
    CTP_LOG(ctp::logDEBUG, *_log)
        << "Using localisation filter for fragment" << fragment
        << " with cutoff " << _loc_threshold << flush;
  }
  if (_use_dQfilter) {
    CTP_LOG(ctp::logDEBUG, *_log)
        << "Using Delta Q filter with cutoff " << _dQ_threshold << flush;
  }
  if (_use_oscfilter && _use_dQfilter) {
    CTP_LOG(ctp::logDEBUG, *_log) << "WARNING: filtering for optically active "
                                     "CT transition - might not make sense... "
                                  << flush;
  }
  if (_use_dQfilter + _use_oscfilter + _use_localisationfilter +
          _use_oscfilter <
      1) {
    CTP_LOG(ctp::logDEBUG, *_log) << "WARNING: No filter is used " << flush;
  }
}

std::vector<int> Statefilter::ComparePairofVectors(
    std::vector<int>& vec1, std::vector<int>& vec2) const {
  std::vector<int> result(std::min(vec1, vec2));
  std::vector<int>::iterator it;
  std::sort(vec1.begin(), vec1.end());
  std::sort(vec2.begin(), vec2.end());
  it = std::set_intersection(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                             result.begin());
  result.resize(it - result.begin());
  return result;
}

std::vector<int> Statefilter::CollapseResults(
    std::vector<std::vector<int> >& results) const {
  if (results.size() == 1) {
    return results[0];
  } else {
    std::vector<int> result = results[0];
    for (unsigned i = 1; i < results.size(); i++) {
      result = ComparePairofVectors(result, results[i]);
    }
    return result;
  }
}

QMState Statefilter::CalcState(Orbitals& orbitals) const {

  if (_use_dQfilter + _use_oscfilter + _use_localisationfilter +
          _use_oscfilter <
      1) {
    return _statehist[0];
  }

  std::vector<std::vector<int> > results;
  if (_use_oscfilter) {
    results.push_back(OscFilter(orbitals));
  }
  if (_use_localisationfilter) {
    results.push_back(LocFilter(orbitals));
  }
  if (_use_overlapfilter) {
    results.push_back(OverlapFilter(orbitals));
  }
  if (_use_dQfilter) {
    results.push_back(DeltaQFilter(orbitals));
  }

  std::vector<int> result = CollapseResults(results);
  QMState state;
  if (result.size() < 1) {
    state = _statehist.back();
    CTP_LOG(ctp::logDEBUG, *_log)
        << "No State found by filter using last state: " << state.ToString()
        << flush;
  } else {
    state = QMState(_statehist.back().Type(), result[0], false);
    CTP_LOG(ctp::logDEBUG, *_log)
        << "Next State is: " << state.ToString() << flush;
  }
  return state;
}

QMState Statefilter::CalcStateAndUpdate(Orbitals& orbitals) {
  QMState result = CalcState(orbitals);
  _statehist.push_back(result);
  if (_use_overlapfilter) {
    UpdateLastCoeff(orbitals);
  }
  return result;
}

std::vector<int> Statefilter::OscFilter(const Orbitals& orbitals) const {
  const std::vector<double> oscs = orbitals.Oscillatorstrengths();
  std::vector<int> indexes;
  for (unsigned i = 0; i < oscs.size(); i++) {
    if (oscs[i] > _oscthreshold) indexes.push_back(i);
  }
  return indexes;
}

std::vector<int> Statefilter::LocFilter(const Orbitals& orbitals) const {
  std::vector<int> indexes;
  const std::vector<Eigen::VectorXd>& popE =
      (_statehist[0].Type() == QMStateType::Singlet)
          ? orbitals.getFragment_E_localisation_singlet()
          : orbitals.getFragment_E_localisation_triplet();
  const std::vector<Eigen::VectorXd>& popH =
      (_statehist[0].Type() == QMStateType::Singlet)
          ? orbitals.getFragment_H_localisation_singlet()
          : orbitals.getFragment_H_localisation_triplet();
  int fragmentindex = 1;
  if (_localiseonA) {
    fragmentindex = 0;
  }
  for (unsigned i = 0; i < popE.size(); i++) {
    if (popE[i](fragmentindex) > _loc_threshold &&
        popH[i](fragmentindex) > _loc_threshold) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

std::vector<int> Statefilter::DeltaQFilter(const Orbitals& orbitals) const {
  std::vector<int> indexes;
  const std::vector<Eigen::VectorXd>& dQ_frag =
      (_statehist[0].Type() == QMStateType::Singlet)
          ? orbitals.getFragmentChargesSingEXC()
          : orbitals.getFragmentChargesTripEXC();
  for (unsigned i = 0; i < dQ_frag.size(); i++) {
    if (std::abs(dQ_frag[i](0)) > _dQ_threshold) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

Eigen::VectorXd Statefilter::CalculateOverlap(Orbitals& orbitals) const {
  BasisSet dftbs;
  dftbs.LoadBasisSet(orbitals.getDFTbasisName());
  AOBasis dftbasis;
  dftbasis.AOBasisFill(dftbs, orbitals.QMAtoms());
  AOOverlap dftoverlap;
  dftoverlap.Fill(dftbasis);
  _S_onehalf = dftoverlap.Pseudo_InvSqrt(1e-8);
  Eigen::MatrixXd ortho_coeffs = CalcOrthoCoeffs(orbitals);
  Eigen::VectorXd overlap = (ortho_coeffs * _laststatecoeff).cwiseAbs2();
  return overlap;
}

Eigen::MatrixXd Statefilter::CalcOrthoCoeffs(Orbitals& orbitals) const {
  QMStateType type = _statehist[0].Type();
  Eigen::MatrixXd ortho_coeffs;
  if (type.isSingleParticleState()) {
    Eigen::MatrixXd coeffs;
    if (type == QMStateType::DQPstate) {
      coeffs = orbitals.CalculateQParticleAORepresentation();
    } else {
      coeffs = orbitals.MOCoefficients();
    }
    ortho_coeffs = _S_onehalf * coeffs;
  } else {
    throw std::runtime_error("Overlap for excitons not implemented yet");
  }
  return ortho_coeffs;
}

void Statefilter::UpdateLastCoeff(Orbitals& orbitals) {
  Eigen::MatrixXd ortho_coeffs = CalcOrthoCoeffs(orbitals);
  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }
  _laststatecoeff = ortho_coeffs.col(_statehist.back().Index() - offset);
}

std::vector<int> Statefilter::OverlapFilter(Orbitals& orbitals) const {
  std::vector<int> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<int>{_statehist[0].Index()};
    return indexes;
  }

  Eigen::VectorXd Overlap = CalculateOverlap(orbitals);
  int validelements = Overlap.size();
  for (int i = 0; i < Overlap.size(); i++) {
    if (Overlap(i) < _overlapthreshold) {
      validelements--;
    }
  }

  std::vector<int> index = std::vector<int>(Overlap.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&Overlap](int i1, int i2) {
    return Overlap[i1] > Overlap[i2];
  });

  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (int i : index) {
    if (int(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}

}  // namespace xtp
}  // namespace votca
