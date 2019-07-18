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

#include <boost/format.hpp>
#include <votca/xtp/esp2multipole.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/nbo.h>
#include <votca/xtp/populationanalysis.h>

namespace votca {
namespace xtp {
using std::flush;

void Esp2multipole::Initialize(tools::Property& options) {
  std::string key = Identify();
  _do_svd = false;

  _use_mulliken = false;
  _use_CHELPG = false;
  _use_lowdin = false;

  std::string statestring =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".state");
  _state.FromString(statestring);

  std::vector<std::string> choices = {"mulliken", "loewdin", "CHELPG"};
  _method = options.ifExistsAndinListReturnElseThrowRuntimeError(
      key + ".method", choices);

  if (_method == "mulliken")
    _use_mulliken = true;
  else if (_method == "loewdin")
    _use_lowdin = true;
  else if (_method == "CHELPG")
    _use_CHELPG = true;

  if (options.exists(key + ".constraints")) {
    if (options.exists(key + ".constraints.regions")) {
      std::vector<tools::Property*> prop_region =
          options.Select(key + ".constraints.regions.region");
      for (tools::Property* prop : prop_region) {
        std::string indices = prop->get("indices").as<std::string>();
        QMFragment<double> reg = QMFragment<double>("Constraint", 0, indices);
        reg.value() = prop->get("charge").as<double>();
        _regionconstraint.push_back(reg);
        XTP_LOG(logDEBUG, _log) << "Fit constrained by Region" << flush;
        XTP_LOG(logDEBUG, _log) << reg;
      }
    }
    if (options.exists(key + ".constraints.pairs")) {
      std::vector<tools::Property*> prop_pair =
          options.Select(key + ".constraints.pairs.pair");
      for (tools::Property* prop : prop_pair) {
        std::string pairstring = prop->as<std::string>();
        tools::Tokenizer tok(pairstring, "\n\t ,");
        std::vector<int> pairvec;
        tok.ConvertToVector<int>(pairvec);
        std::pair<int, int> pair;
        pair.first = pairvec[0];
        pair.second = pairvec[1];
        _pairconstraint.push_back(pair);
        XTP_LOG(logDEBUG, _log) << "Charge " << pair.first << " " << pair.second
                                << " constrained to be equal." << flush;
      }
    }
  }

  _gridsize = options.ifExistsReturnElseReturnDefault<std::string>(
      key + ".gridsize", "medium");
  _openmp_threads =
      options.ifExistsReturnElseReturnDefault<int>(key + ".openmp", 1);

  if (options.exists(key + ".svd")) {
    _do_svd = options.get(key + ".svd.do_svd").as<bool>();
    _conditionnumber = options.get(key + ".svd.conditionnumber").as<double>();
  }

  return;
}

void Esp2multipole::WritetoFile(std::string output_file,
                                const Orbitals& orbitals) {

  std::string data_format = boost::filesystem::extension(output_file);
  if (!(data_format == ".mps")) {
    throw std::runtime_error(
        "Outputfile format not recognized. Export only to .mps");
  }
  std::string tag = "TOOL:" + Identify() + "_" + _state.ToString();

  orbitals.Multipoles().WriteMPS(output_file, tag);
  return;
}

void Esp2multipole::PrintDipoles(Orbitals& orbitals) {
  Eigen::Vector3d classical_dip = orbitals.Multipoles().CalcDipole();

  XTP_LOG(logDEBUG, _log)
      << "El Dipole from fitted charges [e*bohr]:\n\t\t"
      << boost::format(
             " dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f") %
             classical_dip[0] % classical_dip[1] % classical_dip[2] %
             classical_dip.squaredNorm()
      << flush;
  Eigen::Vector3d qm_dip = orbitals.CalcElDipole(_state);
  XTP_LOG(logDEBUG, _log)
      << "El Dipole from exact qm density [e*bohr]:\n\t\t"
      << boost::format(
             " dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f") %
             qm_dip[0] % qm_dip[1] % qm_dip[2] % qm_dip.squaredNorm()
      << flush;
}

void Esp2multipole::Extractingcharges(Orbitals& orbitals) {
  XTP_LOG(logDEBUG, _log) << "===== Running on " << OPENMP::getMaxThreads()
                          << " threads ===== " << flush;

  if (_use_mulliken) {
    Mulliken mulliken;
    mulliken.CalcChargeperAtom(orbitals, _state);
  } else if (_use_lowdin) {
    Lowdin lowdin;
    lowdin.CalcChargeperAtom(orbitals, _state);
  } else if (_use_CHELPG) {
    Espfit esp = Espfit(_log);
    if (_pairconstraint.size() > 0) {
      esp.setPairConstraint(_pairconstraint);
    }
    if (_regionconstraint.size() > 0) {
      esp.setRegionConstraint(_regionconstraint);
    }

    if (_do_svd) {
      esp.setUseSVD(_conditionnumber);
    }
    esp.Fit2Density(orbitals, _state, _gridsize);
  }

  PrintDipoles(orbitals);
}

}  // namespace xtp
}  // namespace votca
