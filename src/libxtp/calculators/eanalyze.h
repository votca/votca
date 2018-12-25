/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_EANALYZE_H
#define _VOTCA_XTP_EANALYZE_H

#include <votca/xtp/qmcalculator.h>
#include <math.h>
#include <votca/tools/tokenizer.h>
#include <votca/xtp/qmstate.h>
#include <votca/tools/histogramnew.h>
#include <fstream>
#include <numeric>

namespace votca {
    namespace xtp {

        class EAnalyze : public QMCalculator {
        public:

            EAnalyze() {
            };

            ~EAnalyze() {
            };

            std::string Identify() {
                return "eanalyze";
            }

            void Initialize(tools::Property *opt);
            bool EvaluateFrame(Topology *top);
            void SiteHist(Topology *top, QMStateType state);
            void PairHist(Topology *top, QMStateType state);
            void SiteCorr(Topology *top, QMStateType state);

        private:

            double _resolution_pairs;
            double _resolution_sites;
            double _resolution_space;
            std::string _distancemode;

            std::vector<QMStateType> _states;

            bool _skip_corr;
            bool _skip_sites;
            bool _skip_pairs;

            bool _doenergy_landscape;
            int _first_seg;
            int _last_seg;

            std::string _seg_pattern;
            std::vector<Segment*> _seg_shortlist;

        };

        void EAnalyze::Initialize(tools::Property *opt) {
            _skip_corr = false;
            _skip_sites = false;
            _skip_pairs = false;
            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults(opt, "xtp");
            std::string key = "options." + Identify();
            if (opt->exists(key + ".resolution_pairs")) {
                _resolution_pairs = opt->get(key + ".resolution_pairs").as< double >();
            } else {
                _skip_pairs = true;
            }
            if (opt->exists(key + ".resolution_sites")) {
                _resolution_sites = opt->get(key + ".resolution_sites").as< double >();
            } else {
                _skip_sites = true;
            }
            if (opt->exists(key + ".resolution_space")) {
                _resolution_space = opt->get(key + ".resolution_space").as< double >();
            } else {
                _skip_corr = true;
            }

            if (opt->exists(key + ".pattern")) {
                _seg_pattern = opt->get(key + ".pattern").as<std::string>();
            } else _seg_pattern = "*";

            if (opt->exists(key + ".states")) {
                std::string statestrings = opt->get(key + ".states").as< std::string>();
                tools::Tokenizer tok(statestrings, ",\n\t ");
                std::vector<std::string> string_vec;
                tok.ToVector(string_vec);
                for (std::string& state : string_vec) {
                    _states.push_back(QMStateType(state));
                }
            } else {
                _states.push_back(QMStateType(QMStateType::Electron));
                _states.push_back(QMStateType(QMStateType::Hole));
            }

            _doenergy_landscape = opt->ifExistsReturnElseReturnDefault<bool>(key + ".energy_landscape", false);


            if (opt->exists(key + ".distancemode")) {
                std::vector<std::string> choices = {"segment", "centerofmass"};
                _distancemode = opt->ifExistsAndinListReturnElseThrowRuntimeError<std::string>(key + ".distancemode", choices);
            } else {
                _distancemode = "segment";
            }

        }

        bool EAnalyze::EvaluateFrame(Topology *top) {

            // Short-list segments according to pattern
            for (Segment* seg : top->Segments()) {
                std::string seg_name = seg->getName();
                if (votca::tools::wildcmp(_seg_pattern.c_str(), seg_name.c_str())) {
                    _seg_shortlist.push_back(seg);
                }
            }
            std::cout << std::endl << "... ... Short-listed " << _seg_shortlist.size()
                    << " segments (pattern='" << _seg_pattern << "')" << std::flush;
            std::cout << std::endl << "... ... ... NOTE Statistics of site energies and spatial"
                    << " correlations thereof are based on the short-listed segments only. "
                    << std::flush;
            std::cout << std::endl << "... ... ...      "
                    << "Statistics of site-energy differences operate on the full list."
                    << std::flush;

            // Calculate
            // ... Site-energy histogram, mean, width
            // ... Pair-energy histogram, mean, width
            // ... Site-energy correlation

            QMNBList &nblist = top->NBList();

            for (QMStateType state : _states) {
                std::cout << std::endl << "... ... excited state " << state.ToString() << std::flush;

                if (!_seg_shortlist.size()) {
                    std::cout << std::endl << "... ... ... No segments short-listed. Skip ... "
                            << std::flush;
                } else {
                    if (_skip_sites) {
                        std::cout << std::endl << "... ... ... Skip site-energy hist." << std::flush;
                    } else {
                        SiteHist(top, state);
                    }
                    if (_skip_corr) {
                        std::cout << std::endl << "... ... ... Skip correlation ..." << std::flush;
                    } else {
                        SiteCorr(top, state);
                    }
                }

                if (!nblist.size()) {
                    std::cout << std::endl << "... ... ... No pairs in topology. Skip ... "
                            << std::flush;
                } else {
                    if (_skip_pairs) {
                        std::cout << std::endl << "... ... ... Skip pair-energy hist." << std::flush;
                    } else {
                        PairHist(top, state);
                    }
                }
            }

            return true;
        }

        void EAnalyze::SiteHist(Topology *top, QMStateType state) {

            std::vector< double > Es;
            Es.reserve(_seg_shortlist.size());
            for (Segment* seg : _seg_shortlist) {
                double E = seg->getSiteEnergy(state.ToSegIndex());
                Es.push_back(E);
            }

            double MAX = *std::max_element(Es.begin(), Es.end());
            double MIN = *std::min_element(Es.begin(), Es.end());
            double sum = std::accumulate(Es.begin(), Es.end(), 0.0);
            double AVG = sum / Es.size();
            double sq_sum = std::inner_product(Es.begin(), Es.end(), Es.begin(), 0.0);
            double STD = std::sqrt(sq_sum / Es.size() - AVG * AVG);

            // Prepare bins
            int BIN = int( (MAX - MIN) / _resolution_sites + 0.5) + 1;

            tools::HistogramNew hist;
            hist.Initialize(MIN, MAX, BIN);
            hist.ProcessRange<std::vector<double>::iterator>(Es.begin(), Es.end());
            tools::Table& tab = hist.data();
            tab.flags()= std::vector<char>(tab.size(), ' ');
            std::string comment = (boost::format("EANALYZE: SITE-ENERGY HISTOGRAM \n # AVG %1$4.7f STD %2$4.7f MIN %3$4.7f MAX %4$4.7f") % AVG % STD % MIN % MAX).str();
            std::string filename = "eanalyze.sitehist_" + state.ToString() + ".out";
            tab.set_comment(comment);
            tab.Save(filename);

            // Write "seg x y z energy" with atomic {x,y,z}
            if (_doenergy_landscape) {
                std::string filename = "eanalyze.landscape_" + state.ToString() + ".out";
                std::ofstream out;
                out.open(filename.c_str());
                if (!out) throw std::runtime_error("error, cannot open file " + filename);
                for (Segment* seg : _seg_shortlist) {
                    if (seg->getId() < _first_seg) {
                        continue;
                    }
                    if (seg->getId() == _last_seg) {
                        break;
                    }
                    double E = seg->getSiteEnergy(state.ToSegIndex());
                    for (Atom *atm : seg->Atoms()) {
                        out << boost::format("%1$3s %2$4.7f %3$4.7f %4$4.7f %5$4.7f\n")
                                % seg->getName()
                                % atm->getPos().getX() % atm->getPos().getY() % atm->getPos().getZ()
                                % E;
                    }
                }
                out.close();
            }
        }

        void EAnalyze::PairHist(Topology *top, QMStateType state) {

            QMNBList &nblist = top->NBList();

            std::string filenamelist = "eanalyze.pairlist_" + state.ToString() + ".out";

            // Collect site-energy differences from neighbourlist
            std::vector< double > dE;
            dE.reserve(2 * nblist.size());
            std::ofstream out;
            out.open(filenamelist.c_str());
            if (!out) throw std::runtime_error("error, cannot open file " + filenamelist);
            for (QMPair *pair : nblist) {
                Segment *seg1 = pair->Seg1();
                Segment *seg2 = pair->Seg2();
                double deltaE = seg2->getSiteEnergy(state.ToSegIndex()) - seg1->getSiteEnergy(state.ToSegIndex());
                dE.push_back(deltaE);
                dE.push_back(-deltaE);
                out << boost::format("%1$5d %2$5d %3$4.7f \n") % seg1->getId() % seg2->getId() % deltaE;
            }
            out.close();

            double MAX = *std::max_element(dE.begin(), dE.end());
            double MIN = *std::min_element(dE.begin(), dE.end());
            double sum = std::accumulate(dE.begin(), dE.end(), 0.0);
            double AVG = sum / dE.size();
            double sq_sum = std::inner_product(dE.begin(), dE.end(), dE.begin(), 0.0);
            double STD = std::sqrt(sq_sum / dE.size() - AVG * AVG);
            int BIN = int( (MAX - MIN) / _resolution_pairs + 0.5) + 1;

            std::string filename2 = "eanalyze.pairhist_" + state.ToString() + ".out";
            tools::HistogramNew hist;
            hist.Initialize(MIN, MAX, BIN);
            hist.ProcessRange<std::vector<double>::iterator>(dE.begin(), dE.end());
            tools::Table& tab = hist.data();
            std::string comment = (boost::format("EANALYZE: PAIR-ENERGY HISTOGRAM \n # AVG %1$4.7f STD %2$4.7f MIN %3$4.7f MAX %4$4.7f") % AVG % STD % MIN % MAX).str();
            tab.set_comment(comment);
            tab.flags()= std::vector<char>(tab.size(), ' ');
            tab.Save(filename2);
        }

        void EAnalyze::SiteCorr(Topology *top, QMStateType state) {

            std::vector< double > Es;
            Es.reserve(_seg_shortlist.size());
            for (Segment* seg : _seg_shortlist) {
                double E = seg->getSiteEnergy(state.ToSegIndex());
                Es.push_back(E);
            }

            double sum = std::accumulate(Es.begin(), Es.end(), 0.0);
            double AVG = sum / Es.size();
            double sq_sum = std::inner_product(Es.begin(), Es.end(), Es.begin(), 0.0);
            double VAR = sq_sum / Es.size() - AVG * AVG;
            double STD = std::sqrt(VAR);

            // Collect inter-site distances, correlation product
            std::vector< Segment* > ::iterator sit1;
            std::vector< Segment* > ::iterator sit2;

            tools::Table tabcorr;
            int length = _seg_shortlist.size()*(_seg_shortlist.size() - 1) / 2;
            tabcorr.resize(length);
            int index = 0;
            for (sit1 = _seg_shortlist.begin(); sit1 < _seg_shortlist.end(); ++sit1) {
                for (sit2 = sit1 + 1; sit2 < _seg_shortlist.end(); ++sit2) {
                    double R = abs(top->PbShortestConnect((*sit1)->getPos(),
                            (*sit2)->getPos()));
                    if (_distancemode == "segment") {
                        for (Fragment* frag1 : (*sit1)->Fragments()) {
                            for (Fragment* frag2 : (*sit2)->Fragments()) {
                                double R_FF = tools::abs(top->PbShortestConnect(frag1->getPos(),
                                        frag2->getPos()));
                                if (R_FF < R) {
                                    R = R_FF;
                                }
                            }
                        }
                    }
                    double C = ((*sit1)->getSiteEnergy(state.ToSegIndex()) - AVG)
                            * ((*sit2)->getSiteEnergy(state.ToSegIndex()) - AVG);
                    tabcorr.set(index, R, C);
                    index++;
                }
            }

            double MIN = tabcorr.x().minCoeff();
            double MAX = tabcorr.x().maxCoeff();
            std::string corrfile = "eanalyze.sitecorr.atomic_" + state.ToString() + ".out";
            tabcorr.Save(corrfile);

            // Prepare bins
            int BIN = int( (MAX - MIN) / _resolution_space + 0.5) + 1;
            std::vector< std::vector<double> > histCs;
            histCs.resize(BIN);

            for (int i = 0; i < tabcorr.size(); ++i) {
                int bin = int((tabcorr.x()[i] - MIN) / _resolution_space + 0.5);
                histCs[bin].push_back(tabcorr.y()[i]);
            }

            tools::Table histC;
            histC.SetHasYErr(true);
            // Calculate spatial correlation
            histC.resize(BIN);
            for (int bin = 0; bin < BIN; ++bin) {
                double corr = 0.0;
                double dcorr2 = 0.0;
                for (unsigned i = 0; i < histCs[bin].size(); ++i) {
                    corr += histCs[bin][i] / VAR;
                }
                corr = corr / histCs[bin].size();
                for (unsigned i = 0; i < histCs[bin].size(); ++i) {
                    dcorr2 += (histCs[bin][i] / VAR / histCs[bin].size() - corr)*(histCs[bin][i] / VAR / histCs[bin].size() - corr);
                }
                // error on mean value
                dcorr2 = dcorr2 / histCs[bin].size() / (histCs[bin].size() - 1);
                double R = MIN + bin*_resolution_space;
                histC.set(bin, R, corr, ' ', std::sqrt(dcorr2));
            }

            std::string filename = "eanalyze.sitecorr_" + state.ToString() + ".out";
            std::string comment = (boost::format("EANALYZE:  SPATIAL SITE-ENERGY CORRELATION \n # AVG %1$4.7f STD %2$4.7f MIN %3$4.7f MAX %4$4.7f") % AVG % STD % MIN % MAX).str();
            histC.set_comment(comment);
            histC.Save(filename);

        }

    }
}

#endif // _VOTCA_XTP_EANALYZE_H
