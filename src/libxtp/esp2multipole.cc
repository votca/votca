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


#include <votca/xtp/esp2multipole.h>
#include <boost/format.hpp>
#include <votca/xtp/orbitals.h>

namespace votca {
    namespace xtp {

        void Esp2multipole::Initialize(Property* options) {
            string key = Identify();
            _use_ecp = false;
            _do_svd = false;

            _use_mulliken = false;
            _use_CHELPG = false;
            _use_bulkESP = false;
            _use_GDMA = false;
            _use_CHELPG_SVD = false;
            _use_lowdin = false;
            _use_NBO = false;

            _state = options->get(key + ".state").as<string> ();
            _state_no = options->get(key + ".statenumber").as<int> ();
            _spin = options->get(key + ".spin").as<string> ();
            if (options->exists(key + ".ecp")) {
                _use_ecp = options->get(key + ".ecp").as<bool> ();
            }
            if (options->exists(key + ".method")) {
                _method = options->get(key + ".method").as<string> ();
                if (_method == "Mulliken" || _method == "mulliken")_use_mulliken = true;
                else if (_method == "loewdin" || _method == "Loewdin") _use_lowdin = true;
                else if (_method == "CHELPG")_use_CHELPG = true;
                else if (_method == "GDMA") throw std::runtime_error("GDMA not implemented yet");
                else if (_method == "CHELPG_SVD") throw std::runtime_error("CHELPG_SVD not implemented yet");
                else if (_method == "NBO") _use_NBO = true;
                else throw std::runtime_error("Method not recognized. Mulliken, Lowdin and CHELPG implemented");
            } else _use_CHELPG = true;

            if (_use_CHELPG) {
                _integrationmethod = options->get(key + ".integrationmethod").as<string> ();
            }
            
            if (options->exists(key + ".constraints")) {
                 if (options->exists(key + ".constraints.regions")) {
                     std::list<Property*> prop_region = options->Select(key + ".constraints.regions.region");
                     for (std::list<Property*> ::iterator it = prop_region.begin(); it != prop_region.end(); ++it) {
                         std::string indices=(*it)->get("indices").as<std::string>();
                         tools::Tokenizer tok(indices,"\n\t ,");
                         Espfit::region reg;
                         tok.ConvertToVector<int>(reg.atomindices);
                         reg.charge=(*it)->get("charge").as<double>();
                         _regionconstraint.push_back(reg);
                         CTP_LOG(ctp::logDEBUG, *_log) << "Fit constrained by SUM(";
                         for(int i:reg.atomindices){
                             CTP_LOG(ctp::logDEBUG, *_log)<<i<<" ";
                         }
                        CTP_LOG(ctp::logDEBUG, *_log)<<")="<<reg.charge<< flush;
                     }
                 }
                 if (options->exists(key + ".constraints.pairs")) {
                     std::list<Property*> prop_pair = options->Select(key + ".constraints.pairs.pair");
                     for (std::list<Property*> ::iterator it = prop_pair.begin(); it != prop_pair.end(); ++it) {
                         std::string pairstring=(*it)->as<std::string>();
                        tools::Tokenizer tok(pairstring,"\n\t ,");
                        std::vector<int> pairvec;
                        tok.ConvertToVector<int>(pairvec);
                        std::pair<int,int> pair;
                        pair.first=pairvec[0];
                        pair.second=pairvec[1];
                        _pairconstraint.push_back(pair);
                        CTP_LOG(ctp::logDEBUG, *_log) << "Charge "<<pair.first<<" "<<pair.second<<" constrained to be equal."<<flush;
                     }
                 }
            }


            if (!(_integrationmethod == "numeric" || _integrationmethod == "analytic")) {
                std::runtime_error("Method not recognized. Only numeric and analytic available");
            }
            if (options->exists(key + ".gridsize")) {
                _gridsize = options->get(key + ".gridsize").as<string>();
            } else _gridsize = "medium";
            if (options->exists(key + ".openmp")) {
                _openmp_threads = options->get(key + ".openmp").as<int>();
            } else _openmp_threads = 0;
            if (options->exists(key + ".svd")) {
                _do_svd = options->get(key + ".svd.do_svd").as<bool>();
                _conditionnumber = options->get(key + ".svd.conditionnumber").as<double>();
            }

            // get the path to the shared folders with xml files
            char *votca_share = getenv("VOTCASHARE");
            if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
            return;
        }

        string Esp2multipole::GetIdentifier() {
            string identifier;
            if (_state == "transition" && _spin == "singlet") {
                identifier = (boost::format("n2s%i") % _state_no).str();
            }
            else if (_state == "transition" && _spin == "triplet") {
                identifier = (boost::format("n2t%i") % _state_no).str();
            }
            else if (_state == "excited" && _spin == "triplet") {
                identifier = (boost::format("t%i") % _state_no).str();
            }
            else if (_state == "excited" && _spin == "singlet") {
                identifier = (boost::format("s%i") % _state_no).str();
            }
            else if (_state == "excited" && _spin == "singlet") {
                identifier = (boost::format("s%i") % _state_no).str();
            }
            else if (_state == "ground") {
                identifier = "n";
            }
            else {
                throw std::runtime_error("Esp2multipole::GetIdentifier did not recognize config");
            }
            return identifier;
        }

        void Esp2multipole::WritetoFile(string _output_file, string identifier) {

            string data_format = boost::filesystem::extension(_output_file);
            if (!(data_format == ".mps")) {
                throw std::runtime_error("Outputfile format not recognized. Export only to .mps");
            }
            string tag = "TOOL:" + Identify() + "_" + GetIdentifier() + "_" + _spin;

            QMInterface Converter;
            ctp::PolarSeg result = Converter.Convert(_Atomlist);

            result.WriteMPS(_output_file, tag);
            return;
        }

        void Esp2multipole::Extractingcharges(Orbitals & _orbitals) {
            int threads = 1;
#ifdef _OPENMP
            if (_openmp_threads > 0) omp_set_num_threads(_openmp_threads);
            threads = omp_get_max_threads();
#endif
            CTP_LOG(ctp::logDEBUG, *_log) << "===== Running on " << threads << " threads ===== " << flush;

            _Atomlist = _orbitals.QMAtoms();
            Eigen::MatrixXd DMAT_tot;
            BasisSet bs;
            bs.LoadBasisSet(_orbitals.getDFTbasis());
            AOBasis basis;
            basis.AOBasisFill(bs, _Atomlist);

            bool _do_transition = false;
            if (_state == "transition") {
              
                _do_transition = true;
                DMAT_tot = _orbitals.TransitionDensityMatrix(_spin, _state_no - 1);

            } else if (_state == "ground" || _state == "excited") {
                Eigen::MatrixXd DMATGS = _orbitals.DensityMatrixGroundState();
                DMAT_tot = DMATGS;
                if (_state_no > 0 && _state == "excited") {
                    std::vector<Eigen::MatrixXd > DMAT = _orbitals.DensityMatrixExcitedState(_spin, _state_no - 1);
                    DMAT_tot = DMAT_tot - DMAT[0] + DMAT[1];
                }
                // Ground state + hole_contribution + electron contribution
            } else throw std::runtime_error("State entry not recognized");

            if (_use_mulliken) {
                Mulliken mulliken;
                mulliken.EvaluateMulliken(_Atomlist, DMAT_tot, basis, _do_transition);
            }
            else if (_use_lowdin) {
                Lowdin lowdin;
                lowdin.EvaluateLowdin(_Atomlist, DMAT_tot, basis, _do_transition);
            } else if (_use_CHELPG) {
                Espfit esp = Espfit(_log);
                if(_pairconstraint.size()>0){
                    esp.setPairConstraint(_pairconstraint);
                }
                if(_regionconstraint.size()>0){
                    esp.setRegionConstraint(_regionconstraint);
                }
                
                if (_do_svd) {
                    esp.setUseSVD(_conditionnumber);
                }
                if (_integrationmethod == "numeric") {
                    esp.Fit2Density(_Atomlist, DMAT_tot, basis, _gridsize);
                } else if (_integrationmethod == "analytic") esp.Fit2Density_analytic(_Atomlist, DMAT_tot, basis);
            } else if (_use_NBO) {
                std::cout << "WARNING: NBO analysis isn't fully implemented yet." << std::endl;
                //CTP_LOG(logDEBUG, _log) << "Initializing NBO" << flush;
                NBO nbo = NBO(_log);
                nbo.EvaluateNBO(_Atomlist, DMAT_tot, basis, bs);
            } else {
                std::cout << "Method not recognized." << std::endl;
            }
        }

    }
}
