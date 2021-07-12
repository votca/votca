/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Local private VOTCA includes
#include "eanalyze.h"
#include "votca/xtp/qmstate.h"

namespace votca {
namespace xtp {

void EAnalyze::ParseOptions(const tools::Property &options) {

  resolution_pairs_ = options.get(".resolution_pairs").as<double>();
  resolution_sites_ = options.get(".resolution_sites").as<double>();
  resolution_spatial_ = options.get(".resolution_spatial").as<double>();
  seg_pattern_ = options.get(".match_pattern").as<std::string>();

  states_ = options.get(".states").as<std::vector<QMStateType>>();

  doenergy_landscape_ = options.get(".output_energy_landscape").as<bool>();

  if (options.exists(".distancemode")) {
    std::string distancemode = options.get("distancemode").as<std::string>();
    if (distancemode == "centerofmass") {
      atomdistances_ = false;
    } else {
      atomdistances_ = true;
    }
  }
}

bool EAnalyze::Evaluate(Topology &top) {

  // Short-list segments according to pattern
  for (Segment &seg : top.Segments()) {
    std::string seg_name = seg.getType();
    if (votca::tools::wildcmp(seg_pattern_, seg_name)) {
      seg_shortlist_.push_back(&seg);
    }
  }
  std::cout << std::endl
            << "... ... Short-listed " << seg_shortlist_.size()
            << " segments (pattern='" << seg_pattern_ << "')" << std::flush;
  std::cout
      << std::endl
      << "... ... ... NOTE Statistics of site energies and spatial"
      << " correlations thereof are based on the short-listed segments only. "
      << std::flush;
  std::cout << std::endl
            << "... ... ...      "
            << "Statistics of site-energy differences operate on the full list."
            << std::flush;

  // Calculate
  // ... Site-energy histogram, mean, width
  // ... Pair-energy histogram, mean, width
  // ... Site-energy correlation

  const QMNBList &nblist = top.NBList();

  for (QMStateType state : states_) {
    std::cout << std::endl
              << "... ... excited state " << state.ToString() << std::flush;

    if (!seg_shortlist_.size()) {
      std::cout << std::endl
                << "... ... ... No segments short-listed. Skip ... "
                << std::flush;
    } else {
      SiteHist(state);
      SiteCorr(top, state);
    }

    if (!nblist.size()) {
      std::cout << std::endl
                << "... ... ... No pairs in topology. Skip ... " << std::flush;
    } else {
      PairHist(top, state);
    }
  }

  return true;
}

void EAnalyze::SiteHist(QMStateType state) const {

  std::vector<double> Es;
  Es.reserve(seg_shortlist_.size());
  for (Segment *seg : seg_shortlist_) {
    Es.push_back(seg->getSiteEnergy(state) * tools::conv::hrt2ev);
  }

  double MAX = *std::max_element(Es.begin(), Es.end());
  double MIN = *std::min_element(Es.begin(), Es.end());
  double sum = std::accumulate(Es.begin(), Es.end(), 0.0);
  double AVG = sum / double(Es.size());
  double sq_sum = std::inner_product(Es.begin(), Es.end(), Es.begin(), 0.0);
  double STD = std::sqrt(sq_sum / double(Es.size()) - AVG * AVG);

  // Prepare bins
  Index BIN = Index((MAX - MIN) / resolution_sites_ + 0.5) + 1;

  tools::HistogramNew hist;
  hist.Initialize(MIN, MAX, BIN);
  hist.ProcessRange<std::vector<double>::iterator>(Es.begin(), Es.end());
  tools::Table &tab = hist.data();
  tab.flags() = std::vector<char>(tab.size(), ' ');
  std::string comment =
      (boost::format("EANALYZE: SITE-ENERGY HISTOGRAM[eV] \n # AVG %1$4.7f STD "
                     "%2$4.7f MIN %3$4.7f MAX %4$4.7f") %
       AVG % STD % MIN % MAX)
          .str();
  std::string filename = "eanalyze.sitehist_" + state.ToString() + ".out";
  tab.set_comment(comment);
  tab.Save(filename);

  // Write "seg x y z energy" with atomic {x,y,z}
  if (doenergy_landscape_) {
    std::string filename2 = "eanalyze.landscape_" + state.ToString() + ".out";
    std::ofstream out;
    out.open(filename2);
    if (!out) {
      throw std::runtime_error("error, cannot open file " + filename2);
    }
    for (Segment *seg : seg_shortlist_) {
      if (seg->getId() < first_seg_) {
        continue;
      }
      if (seg->getId() == last_seg_) {
        break;
      }
      double E = seg->getSiteEnergy(state);
      for (Atom &atm : *seg) {
        out << boost::format("%1$3s %2$4.7f %3$4.7f %4$4.7f %5$4.7f\n") %
                   seg->getType() % atm.getPos().x() % atm.getPos().y() %
                   atm.getPos().z() % E;
      }
    }
    out.close();
  }
}

void EAnalyze::PairHist(const Topology &top, QMStateType state) const {

  const QMNBList &nblist = top.NBList();

  std::string filenamelist = "eanalyze.pairlist_" + state.ToString() + ".out";

  // Collect site-energy differences from neighbourlist
  std::vector<double> dE;
  dE.reserve(2 * nblist.size());
  std::ofstream out;
  out.open(filenamelist);
  if (!out) {
    throw std::runtime_error("error, cannot open file " + filenamelist);
  }
  for (QMPair *pair : nblist) {
    double deltaE = pair->getdE12(state) * tools::conv::hrt2ev;
    dE.push_back(deltaE);
    dE.push_back(-deltaE);
    out << boost::format("%1$5d %2$5d %3$4.7f \n") % pair->Seg1()->getId() %
               pair->Seg2()->getId() % deltaE;
  }
  out.close();

  double MAX = *std::max_element(dE.begin(), dE.end());
  double MIN = *std::min_element(dE.begin(), dE.end());
  double sum = std::accumulate(dE.begin(), dE.end(), 0.0);
  double AVG = sum / double(dE.size());
  double sq_sum = std::inner_product(dE.begin(), dE.end(), dE.begin(), 0.0);
  double STD = std::sqrt(sq_sum / double(dE.size()) - AVG * AVG);
  Index BIN = Index((MAX - MIN) / resolution_pairs_ + 0.5) + 1;

  std::string filename2 = "eanalyze.pairhist_" + state.ToString() + ".out";
  tools::HistogramNew hist;
  hist.Initialize(MIN, MAX, BIN);
  hist.ProcessRange<std::vector<double>::iterator>(dE.begin(), dE.end());
  tools::Table &tab = hist.data();
  std::string comment =
      (boost::format("EANALYZE: PAIR-ENERGY HISTOGRAM[eV] \n # AVG %1$4.7f STD "
                     "%2$4.7f MIN %3$4.7f MAX %4$4.7f") %
       AVG % STD % MIN % MAX)
          .str();
  tab.set_comment(comment);
  tab.flags() = std::vector<char>(tab.size(), ' ');
  tab.Save(filename2);
}

void EAnalyze::SiteCorr(const Topology &top, QMStateType state) const {

  std::vector<double> Es;
  Es.reserve(seg_shortlist_.size());
  for (Segment *seg : seg_shortlist_) {
    Es.push_back(seg->getSiteEnergy(state) * tools::conv::hrt2ev);
  }

  double sum = std::accumulate(Es.begin(), Es.end(), 0.0);
  double AVG = sum / double(Es.size());
  double sq_sum = std::inner_product(Es.begin(), Es.end(), Es.begin(), 0.0);
  double VAR = sq_sum / double(Es.size()) - AVG * AVG;
  double STD = std::sqrt(VAR);

  // Collect inter-site distances, correlation product
  tools::Table tabcorr;
  Index length =
      Index(seg_shortlist_.size()) * Index(seg_shortlist_.size() - 1) / 2;
  tabcorr.resize(length);
  Index index = 0;
  for (Index i = 0; i < Index(seg_shortlist_.size()); i++) {
    const Segment &segi = *seg_shortlist_[i];
    for (Index j = i + 1; j < Index(seg_shortlist_.size()); j++) {
      const Segment &segj = *seg_shortlist_[j];
      double R = 0;
      if (atomdistances_) {
        R = top.GetShortestDist(segi, segj);
      } else {
        R = (top.PbShortestConnect(segi.getPos(), segj.getPos())).norm();
      }

      double C =
          (segi.getSiteEnergy(state) - AVG) * (segj.getSiteEnergy(state) - AVG);
      tabcorr.set(index, R * tools::conv::bohr2nm, C);
      index++;
    }
  }

  double MIN = tabcorr.x().minCoeff();
  double MAX = tabcorr.x().maxCoeff();

  // Prepare bins
  Index BIN = Index((MAX - MIN) / resolution_spatial_ + 0.5) + 1;
  std::vector<std::vector<double>> histCs;
  histCs.resize(BIN);

  for (Index i = 0; i < tabcorr.size(); ++i) {
    Index bin = Index((tabcorr.x()[i] - MIN) / resolution_spatial_ + 0.5);
    histCs[bin].push_back(tabcorr.y()[i]);
  }

  tools::Table histC;
  histC.SetHasYErr(true);
  // Calculate spatial correlation
  histC.resize(BIN);
  for (Index bin = 0; bin < BIN; ++bin) {
    double corr = 0.0;
    double dcorr2 = 0.0;
    for (double entry : histCs[bin]) {
      corr += entry / VAR;
    }
    corr = corr / double(histCs[bin].size());
    for (Index i = 0; i < Index(histCs[bin].size()); ++i) {
      dcorr2 += (histCs[bin][i] / VAR / double(histCs[bin].size()) - corr) *
                (histCs[bin][i] / VAR / double(histCs[bin].size()) - corr);
    }
    // error on mean value
    dcorr2 =
        dcorr2 / double(histCs[bin].size()) / double(histCs[bin].size() - 1);
    double R = MIN + double(bin) * resolution_spatial_;
    histC.set(bin, R, corr, ' ', std::sqrt(dcorr2));
  }

  std::string filename = "eanalyze.sitecorr_" + state.ToString() + ".out";
  std::string comment =
      (boost::format("EANALYZE: DISTANCE[nm] SPATIAL SITE-ENERGY "
                     "CORRELATION[eV] \n # AVG "
                     "%1$4.7f STD %2$4.7f MIN %3$4.7f MAX %4$4.7f") %
       AVG % STD % MIN % MAX)
          .str();
  histC.set_comment(comment);
  histC.Save(filename);
}

}  // namespace xtp
}  // namespace votca
