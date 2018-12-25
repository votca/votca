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

#ifndef VOTCA_XTP_IANALYZE_H
#define VOTCA_XTP_IANALYZE_H

#include <votca/xtp/qmcalculator.h>
#include <math.h>
#include <votca/xtp/qmpair.h>
#include <votca/tools/histogramnew.h>
#include <votca/xtp/qmstate.h>
#include <numeric>

namespace votca {
    namespace xtp {

        class IAnalyze : public QMCalculator {
        public:

            std::string Identify() {
                return "ianalyze";
            }

            void Initialize(tools::Property *options);
            bool EvaluateFrame(Topology *top);
            void IHist(Topology *top, QMStateType state);
            void IRdependence(Topology *top, QMStateType state);

        private:

            double _resolution_logJ2;
            std::vector<QMStateType> _states;
            double _resolution_space;
            std::vector<QMPair::PairType> _pairtype;
            bool _do_pairtype;
            bool _do_IRdependence;

        };

        void IAnalyze::Initialize(tools::Property *opt) {
            std::cout<<std::endl;
            _do_pairtype = false;
            _do_IRdependence = false;
            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults(opt, "xtp");
            std::string key = "options." + Identify();
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

            _resolution_logJ2 = opt->get(key + ".resolution_logJ2").as< double >();
            if (opt->exists(key + ".pairtype")) {
                _do_pairtype = true;
                std::string _store_stdstring = opt->get(key + ".pairtype").as<std::string> ();
                if (_store_stdstring.find("Hopping") != std::string::npos) _pairtype.push_back(QMPair::Hopping);
                if (_store_stdstring.find("Excitoncl") != std::string::npos) _pairtype.push_back(QMPair::Excitoncl);
                if (!_pairtype.size()) {
                    std::cout << std::endl << "... ... No pairtypes recognized will output all pairs. ";
                    _do_pairtype = false;
                }
            }
            if (opt->exists(key + ".resolution_space")) {
                _resolution_space = opt->get(key + ".resolution_space").as< double >();
                if (_resolution_space != 0.0) _do_IRdependence = true;
            }
        }

        bool IAnalyze::EvaluateFrame(Topology *top) {
            std::cout<<std::endl;
            QMNBList &nblist = top->NBList();
            if (!nblist.size()) {
                std::cout << std::endl << "... ... No pairs in topology. Skip...";
                return 0;
            }
            if (_do_pairtype) {
                bool pairs_exist = false;
                for (QMPair* pair:nblist) {
                    QMPair::PairType pairtype = pair->getType();
                    if (std::find(_pairtype.begin(), _pairtype.end(), pairtype) != _pairtype.end()) {
                        pairs_exist = true;
                        break;
                    }
                }
                if (!pairs_exist) {
                    std::cout << std::endl << "... ... No pairs of given pairtypes in topology. Skip...";
                    return 0;
                }
            }
            for (QMStateType state:_states) {
                std::cout<<"Calculating for state "<<state.ToString()<<" now."<<std::endl;
                this->IHist(top, state);
                if (_do_IRdependence) {
                    this->IRdependence(top, state);
                }
            }
            return true;
        }

        void IAnalyze::IHist(Topology *top, QMStateType state) {
            QMNBList &nblist = top->NBList();

            // Collect J2s from pairs
            std::vector< double > J2s;
            for (QMPair* pair:nblist) {
                if (_do_pairtype) {
                    QMPair::PairType pairtype =pair->getType();
                    if (!(std::find(_pairtype.begin(), _pairtype.end(), pairtype) != _pairtype.end())) {
                        continue;
                    }
                }
                double test = pair->getJeff2(state.ToSegIndex());
                if (test <= 0) {
                    continue;
                } // avoid -inf in output
                double J2 = std::log10(test);
                J2s.push_back(J2);
            }

            if (J2s.size() < 1) {
                std::cout<<"WARNING:"+state.ToLongString()+" Couplings are all zero. You have not yet imported them! "<<std::endl;
                return;
            }
            
            double MAX = *std::max_element(J2s.begin(), J2s.end());
            double MIN = *std::min_element(J2s.begin(), J2s.end());
            double sum = std::accumulate(J2s.begin(), J2s.end(), 0.0);
            double AVG = sum / J2s.size();
            double sq_sum = std::inner_product(J2s.begin(), J2s.end(), J2s.begin(), 0.0);
            double STD = std::sqrt(sq_sum / J2s.size() - AVG * AVG);
            // Prepare bins
            int BIN = ((MAX - MIN) / _resolution_logJ2 + 0.5) + 1;
            
            tools::HistogramNew hist;
            hist.Initialize(MIN, MAX, BIN);
            hist.ProcessRange<std::vector<double>::iterator>(J2s.begin(), J2s.end());
            tools::Table& tab = hist.data();
            std::string comment = (boost::format("IANALYZE: PAIR-INTEGRAL J2 HISTOGRAM \n # AVG %1$4.7f STD %2$4.7f MIN %3$4.7f MAX %4$4.7f") % AVG % STD % MIN % MAX).str();
            std::string filename = "ianalyze.ihist_" + state.ToString() + ".out";
            tab.set_comment(comment);
            tab.flags()= std::vector<char>(tab.size(), ' ');
            tab.Save(filename);

        }

        void IAnalyze::IRdependence(Topology *top, QMStateType state) {

            QMNBList &nblist = top->NBList();

            // Collect J2s from pairs
            std::vector< double > J2s;
            J2s.reserve(nblist.size());
            std::vector< double > distances;
            distances.reserve(nblist.size());

            for (QMPair* pair:nblist) {
                double J2 = std::log10(pair->getJeff2(state.ToSegIndex()));
                double distance = tools::abs(pair->R());
                distances.push_back(distance);
                J2s.push_back(J2);
            }
            
            double MAXR = *std::max_element(distances.begin(), distances.end());
            double MINR = *std::min_element(distances.begin(), distances.end());

            // Prepare R bins
            int pointsR = (MAXR - MINR) / _resolution_space;
            std::vector< std::vector<double> > rJ2;
            rJ2.resize(pointsR);

            // Loop over distance
            for (int i = 0; i < pointsR; ++i) {
                double thisMINR = MINR + i*_resolution_space;
                double thisMAXR = MINR + (i + 1) * _resolution_space;
                // now count Js that lie within this R range
                for (unsigned j=0;j<J2s.size();++j) {
                    if (thisMINR < distances[j] && distances[j] < thisMAXR) {
                        rJ2[i].push_back(J2s[j]);
                    }
                }
            }

            tools::Table tab;
            tab.SetHasYErr(true);
            tab.resize(pointsR);
            
            // make plot values
            for (unsigned i=0;i<rJ2.size();i++) {
                const std::vector<double>& vec=rJ2[i];
                double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
                double AVG = sum / vec.size();
                double thisR = MINR + (i + 0.5) * _resolution_space;
                double sq_sum = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
                double STD = std::sqrt(sq_sum / vec.size() - AVG * AVG);
                tab.set(i,thisR,AVG,' ',STD);
            }
           std::string filename = "ianalyze.ispatial_" + state.ToString() + ".out";
           std::string comment ="# IANALYZE: SPATIAL DEPENDENCE OF log10(J2) [r,log10(J),error]";
           tab.setErrorDetails(comment);
           tab.Save(filename);

        }

    }
}



#endif
