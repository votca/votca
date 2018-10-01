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
#include <votca/xtp/qminterface.h>

namespace votca {
    namespace xtp {
      using std::flush;

        void Esp2multipole::Initialize(tools::Property& options) {
            std::string key = Identify();
            _do_svd = false;

            _use_mulliken = false;
            _use_CHELPG = false;
            _use_lowdin = false;

            std::string statestring = options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".state");
            _state.FromString(statestring);
          
            std::vector<std::string> choices={"mulliken","loewdin","CHELPG"};
            _method=options.ifExistsAndinListReturnElseThrowRuntimeError(key+".method",choices);

            if (_method == "mulliken")_use_mulliken = true;
            else if (_method == "loewdin" ) _use_lowdin = true;
            else if (_method == "CHELPG")_use_CHELPG = true;
              

            if (_use_CHELPG) {
                _integrationmethod = options.ifExistsReturnElseReturnDefault<std::string>(key + ".integrationmethod","numeric");
            }
            if (!(_integrationmethod == "numeric" || _integrationmethod == "analytic")) {
                std::runtime_error("Method not recognized. Only numeric and analytic available");
            }
            
            if (options.exists(key + ".constraints")) {
                 if (options.exists(key + ".constraints.regions")) {
                     std::list<tools::Property*> prop_region = options.Select(key + ".constraints.regions.region");
                     for (tools::Property* prop:prop_region) {
                         std::string indices=prop->get("indices").as<std::string>();
                         tools::Tokenizer tok(indices,"\n\t ,");
                         Espfit::region reg;
                         tok.ConvertToVector<int>(reg.atomindices);
                         reg.charge=prop->get("charge").as<double>();
                         _regionconstraint.push_back(reg);
                         XTP_LOG(logDEBUG, *_log) << "Fit constrained by SUM(";
                         for(int i:reg.atomindices){
                             XTP_LOG(logDEBUG, *_log)<<i<<" ";
                         }
                        XTP_LOG(logDEBUG, *_log)<<")="<<reg.charge<< flush;
                     }
                 }
                 if (options.exists(key + ".constraints.pairs")) {
                     std::list<tools::Property*> prop_pair = options.Select(key + ".constraints.pairs.pair");
                     for (tools::Property* prop:prop_pair) {
                         std::string pairstring=prop->as<std::string>();
                        tools::Tokenizer tok(pairstring,"\n\t ,");
                        std::vector<int> pairvec;
                        tok.ConvertToVector<int>(pairvec);
                        std::pair<int,int> pair;
                        pair.first=pairvec[0];
                        pair.second=pairvec[1];
                        _pairconstraint.push_back(pair);
                        XTP_LOG(logDEBUG, *_log) << "Charge "<<pair.first<<" "<<pair.second<<" constrained to be equal."<<flush;
                     }
                 }
            }

          _gridsize =options.ifExistsReturnElseReturnDefault<std::string>(key + ".gridsize","medium");
          _openmp_threads = options.ifExistsReturnElseReturnDefault<int>(key + ".openmp",1);

            if (options.exists(key + ".svd")) {
                _do_svd = options.get(key + ".svd.do_svd").as<bool>();
                _conditionnumber = options.get(key + ".svd.conditionnumber").as<double>();
            }

            // get the path to the shared folders with xml files
            char *votca_share = getenv("VOTCASHARE");
            if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
            return;
        }

        

        void Esp2multipole::WritetoFile(std::string output_file) {

            std::string data_format = boost::filesystem::extension(output_file);
            if (!(data_format == ".mps")) {
                throw std::runtime_error("Outputfile format not recognized. Export only to .mps");
            }
            std::string tag = "TOOL:" + Identify() + "_" + _state.ToString();

            QMInterface Converter;
            PolarSeg result = Converter.Convert(_Atomlist);

            result.WriteMPS(output_file, tag);
            return;
        }

        void Esp2multipole::PrintDipoles(Orbitals& orbitals){
          Eigen::Vector3d CoM = _atomlist.getPos();
          
          Eigen::Vector3d classical_dip = Eigen::Vector3d::Zero();
          for (const QMAtom& atom : _atomlist) {
            classical_dip += (atom.getPos()- CoM) * atom.getPartialcharge();
          }
          CTP_LOG(ctp::logDEBUG, *_log) << "El Dipole from fitted charges [e*bohr]:\n\t\t" 
                  << boost::format(" dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f")
                  % classical_dip[0] % classical_dip[1] % classical_dip[2] % classical_dip.squaredNorm()<< flush;
          Eigen::Vector3d qm_dip=orbitals.CalcElDipole(_state);
          CTP_LOG(ctp::logDEBUG, *_log) << "El Dipole from exact qm density [e*bohr]:\n\t\t"   
                  << boost::format(" dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f")
                  % qm_dip[0] % qm_dip[1] % qm_dip[2] % qm_dip.squaredNorm()<< flush;
        }

        void Esp2multipole::Extractingcharges(Orbitals & orbitals) {
            int threads = 1;
#ifdef _OPENMP
            if (_openmp_threads > 0) omp_set_num_threads(_openmp_threads);
            threads = omp_get_max_threads();
#endif
            XTP_LOG(logDEBUG, *_log) << "===== Running on " << threads << " threads ===== " << flush;

            _atomlist = orbitals.QMAtoms();
            BasisSet bs;
            bs.LoadBasisSet(orbitals.getDFTbasis());
            AOBasis basis;
            basis.AOBasisFill(bs, _atomlist);
            Eigen::MatrixXd DMAT=orbitals.DensityMatrixFull(_state);
            

            if (_use_mulliken) {
                Mulliken mulliken;
                mulliken.EvaluateMulliken(_atomlist, DMAT, basis, _state.isTransition());
            }
            else if (_use_lowdin) {
                Lowdin lowdin;
                lowdin.EvaluateLowdin(_atomlist, DMAT, basis, _state.isTransition());
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
                    esp.Fit2Density(_atomlist, DMAT, basis, _gridsize);
                } else if (_integrationmethod == "analytic") esp.Fit2Density_analytic(_atomlist, DMAT, basis);
            } 

            PrintDipoles(orbitals);
            
        }

    }
}
