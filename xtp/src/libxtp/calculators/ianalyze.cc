/*
 * Copyright 2009-2025 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <cmath>
#include <numeric>

// VOTCA includes
#include <votca/tools/histogram.h>

// Local VOTCA includes
#include "votca/xtp/qmstate.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "ianalyze.h"

namespace votca {
namespace xtp {

void IAnalyze::ParseOptions(const tools::Property &options) {

  states_ = options.get(".states").as<std::vector<QMStateType>>();

  resolution_logJ2_ = options.get(".resolution_logJ2").as<double>();
  if (options.get(".do_pairtype").as<bool>()) {
    do_pairtype_ = true;
    std::string store_stdstring_ = options.get(".pairtype").as<std::string>();
    if (store_stdstring_.find("Hopping") != std::string::npos) {
      pairtype_.push_back(QMPair::Hopping);
    }
    if (store_stdstring_.find("Excitoncl") != std::string::npos) {
      pairtype_.push_back(QMPair::Excitoncl);
    }
    if (!pairtype_.size()) {
      std::cout << "\n... ... No pairtypes recognized will output all pairs. ";
      do_pairtype_ = false;
    }
  }
  if (options.get(".do_resolution_spatial").as<bool>()) {
    resolution_spatial_ = options.get(".resolution_spatial").as<double>();
    if (resolution_spatial_ != 0.0) {
      do_IRdependence_ = true;
    }
  }
}

bool IAnalyze::Evaluate(Topology &top) {
  std::cout << std::endl;
  QMNBList &nblist = top.NBList();
  if (!nblist.size()) {
    std::cout << std::endl << "... ... No pairs in topology. Skip...";
    return false;
  }
  if (do_pairtype_) {
    bool pairs_exist = false;
    for (QMPair *pair : nblist) {
      QMPair::PairType pairtype = pair->getType();
      if (std::find(pairtype_.begin(), pairtype_.end(), pairtype) !=
          pairtype_.end()) {
        pairs_exist = true;
        break;
      }
    }
    if (!pairs_exist) {
      std::cout << std::endl
                << "... ... No pairs of given pairtypes in topology. Skip...";
      return 0;
    }
  }
  for (QMStateType state : states_) {
    std::cout << "Calculating for state " << state.ToString() << " now."
              << std::endl;
    this->IHist(top, state);
    if (do_IRdependence_) {
      IRdependence(top, state);
    }
  }
  return true;
}

void IAnalyze::IHist(Topology &top, QMStateType state) {
  QMNBList &nblist = top.NBList();

  // Collect J2s from pairs
  std::vector<double> J2s;
  for (QMPair *pair : nblist) {
    if (do_pairtype_) {
      QMPair::PairType pairtype = pair->getType();
      if (!(std::find(pairtype_.begin(), pairtype_.end(), pairtype) !=
            pairtype_.end())) {
        continue;
      }
    }
    double test =
        pair->getJeff2(state) * tools::conv::hrt2ev * tools::conv::hrt2ev;
    if (test <= 0) {
      continue;
    }  // avoid -inf in output
    double J2 = std::log10(test);
    J2s.push_back(J2);
  }

  if (J2s.size() < 1) {
    std::cout << "WARNING:" + state.ToLongString() +
                     " Couplings are all zero. You have not yet imported them! "
              << std::endl;
    return;
  }

  double MAX = *std::max_element(J2s.begin(), J2s.end());
  double MIN = *std::min_element(J2s.begin(), J2s.end());
  double sum = std::accumulate(J2s.begin(), J2s.end(), 0.0);
  double AVG = sum / double(J2s.size());
  double sq_sum = std::inner_product(J2s.begin(), J2s.end(), J2s.begin(), 0.0);
  double STD = std::sqrt(sq_sum / double(J2s.size()) - AVG * AVG);
  // Prepare bins
  Index BIN = Index((MAX - MIN) / resolution_logJ2_ + 0.5) + 1;

  tools::Histogram hist;
  hist.Initialize(MIN, MAX, BIN);
  hist.ProcessRange<std::vector<double>::iterator>(J2s.begin(), J2s.end());
  tools::Table &tab = hist.data();
  std::string comment =
      (boost::format("IANALYZE: PAIR-INTEGRAL J2 HISTOGRAM \n # AVG %1$4.7f "
                     "STD %2$4.7f MIN %3$4.7f MAX %4$4.7f") %
       AVG % STD % MIN % MAX)
          .str();
  std::string filename = "ianalyze.ihist_" + state.ToString() + ".out";
  tab.set_comment(comment);
  tab.flags() = std::vector<char>(tab.size(), ' ');
  tab.Save(filename);
}

void IAnalyze::IRdependence(Topology &top, QMStateType state) {

  QMNBList &nblist = top.NBList();

  // Collect J2s from pairs
  std::vector<double> J2s;
  J2s.reserve(nblist.size());
  std::vector<double> distances;
  distances.reserve(nblist.size());

  for (QMPair *pair : nblist) {
    double test =
        pair->getJeff2(state) * tools::conv::hrt2ev * tools::conv::hrt2ev;
    if (test <= 0) {
      continue;
    }  // avoid -inf in output
    double J2 = std::log10(test);
    double distance = pair->R().norm() * tools::conv::bohr2nm;
    distances.push_back(distance);
    J2s.push_back(J2);
  }

  double MAXR = *std::max_element(distances.begin(), distances.end());
  double MINR = *std::min_element(distances.begin(), distances.end());

  // Prepare R bins
  Index pointsR = Index((MAXR - MINR) / resolution_spatial_);
  std::vector<std::vector<double>> rJ2;
  rJ2.resize(pointsR);

  // Loop over distance
  for (Index i = 0; i < pointsR; ++i) {
    double thisMINR = MINR + double(i) * resolution_spatial_;
    double thisMAXR = MINR + double(i + 1) * resolution_spatial_;
    // now count Js that lie within this R range
    for (Index j = 0; j < Index(J2s.size()); ++j) {
      if (thisMINR < distances[j] && distances[j] < thisMAXR) {
        rJ2[i].push_back(J2s[j]);
      }
    }
  }

  tools::Table tab;
  tab.SetHasYErr(true);
  tab.resize(pointsR);

  // make plot values
  for (Index i = 0; i < Index(rJ2.size()); i++) {
    const std::vector<double> &vec = rJ2[i];
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    double AVG = sum / double(vec.size());
    double thisR = MINR + (double(i) + 0.5) * resolution_spatial_;
    double sq_sum =
        std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
    double STD = std::sqrt(sq_sum / double(vec.size()) - AVG * AVG);
    tab.set(i, thisR, AVG, ' ', STD);
  }
  std::string filename = "ianalyze.ispatial_" + state.ToString() + ".out";
  std::string comment =
      "# IANALYZE: SPATIAL DEPENDENCE OF log10(J2) "
      "[r[nm],log10(J)[eV^2],error]";
  tab.setErrorDetails(comment);
  tab.Save(filename);
}

}  // namespace xtp
}  // namespace votca
