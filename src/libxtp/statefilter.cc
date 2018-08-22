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

#include <votca/xtp/statefilter.h>

namespace votca {
    namespace xtp {
      using std::flush;
      
   void Statefilter::Initialize(tools::Property& options){
        if (options.exists("oscillator_strength")) {
          _use_oscfilter = true;
          _oscthreshold = options.ifExistsReturnElseThrowRuntimeError<double>("oscillator_strength");
        }
        if (options.exists("overlap")) {
          _use_overlapfilter = true;
          _overlapthreshold = options.ifExistsReturnElseThrowRuntimeError<double>("overlap");
        }
        if (options.exists("localisation")) {
          _use_localisationfilter = true;
          std::string temp = options.get("localisation").as<std::string> ();
          tools::Tokenizer tok_cleanup(temp, ", \n\t");
          std::vector <std::string> strings_vec;
          tok_cleanup.ToVector(strings_vec);
          if (strings_vec.size() != 2) {
            throw std::runtime_error("qmmmachine: Fragment and localisation threshold are not separated");
          }
          if (strings_vec[0] == "a" || strings_vec[0] == "A") {
            _localiseonA = true;
          } else if (strings_vec[0] == "b" || strings_vec[0] == "B") {
            _localiseonA = false;
          } else {
            throw std::runtime_error("statefiler: Fragment label not known, either A or B");
          }
          _loc_threshold = boost::lexical_cast<double>(strings_vec[1]);
        }

        if (options.exists("charge_transfer")) {
          _use_dQfilter = true;
          _dQ_threshold = options.get("charge_transfer").as<double> ();
        }
        if (_use_dQfilter && _use_localisationfilter) {
          throw std::runtime_error("Cannot use localisation and charge_transfer filter at the same time.");
        }
        if (_use_oscfilter && _use_dQfilter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- WARNING: filtering for optically active CT transition - might not make sense... " << flush;
          }
        
        if(_use_dQfilter+_use_oscfilter+_use_localisationfilter+_oscthreshold<1){
          throw std::runtime_error("No filter is used.");
        }
     
   
   }
   
   void Statefilter::PrintInfo()const{
     CTP_LOG(ctp::logDEBUG, *_log) << " Initial state "<<_statehist[0].ToString() << flush;
     if(_statehist.size()>1){
     CTP_LOG(ctp::logDEBUG, *_log) << " Last state "<<_state.ToString() << flush;
     }
     if(_use_oscfilter){
       CTP_LOG(ctp::logDEBUG, *_log) << " Using oscillator strength filter with cutoff "<<_oscthreshold << flush;
     }
     if(_use_overlapfilter){
       CTP_LOG(ctp::logDEBUG, *_log) << " Using overlap filer with cutoff "<<_overlapthreshold << flush;
     }
     if(_use_localisationfilter){
       std::string fragment="A";
       if(!_localiseonA){
         fragment="B";
       }
       CTP_LOG(ctp::logDEBUG, *_log) << " Using localisation filter for fragment"<<fragment<<" with cutoff "<<_loc_threshold << flush;
     }
     if(_use_dQfilter){
       CTP_LOG(ctp::logDEBUG, *_log) << " Using Delta Q filter with cutoff "<<_dQ_threshold << flush;
     }
     
   }
   
   
   std::vector<int> Statefilter::ComparePairofVectors(std::vector<int> & vec1,std::vector<int> & vec2)const{
    std::vector<int> result(std::min(vec1,vec2));
    std::vector<int>::iterator it;
    std::sort (vec1.begin(), vec1.end());
    std::sort (vec2.begin(), vec2.end());
    it=std::set_intersection (vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), result.begin());
    result.resize(it-result.begin());   
    return result;
   }
   
   std::vector<int> Statefilter::CollapseResults(std::vector< std::vector<int> >& results)const{
     if(results.size()==1){
       return results[0];
     }else{
      std::vector<int> result=results[0];
      for(unsigned i=1;i<results.size();i++){
        result=ComparePairofVectors(result,results[i]);
      }
      return result;
     }
   }
   
   void Statefilter::Filter(const Orbitals& orbitals){
      
     std::vector< std::vector<int> > results;
     if(_use_oscfilter){
       results.push_back(OscFilter(orbitals));
     }
     if(_use_localisationfilter){
       results.push_back(LocFilter(orbitals));
     }
     if(_use_overlapfilter){
       results.push_back(OverlapFilter(orbitals));
     }
     if(_use_dQfilter){
       results.push_back(DeltaQFilter(orbitals));
     }
    _statehist.push_back(_state);
     
    std::vector<int> result=CollapseResults(results);
    if(result.size()<1){
      CTP_LOG(ctp::logDEBUG, *_log) << " No State found by filter using last state"<<_statehist.back().ToString()<< flush;
      _state=_statehist.back();
    }else{
      _state=QMState(_initial_state.Type(),result[0],false);
      CTP_LOG(ctp::logDEBUG, *_log) << " Next State is "<<_state.ToString()<< flush;
    }
   }
   
   std::vector<int> Statefilter::OscFilter(const Orbitals& orbitals){
     const std::vector<double>oscs = orbitals.Oscillatorstrengths();
     std::vector<int> indexes;
     for (unsigned i = 0; i < oscs.size(); i++) {
              if (oscs[i] > _oscthreshold) indexes.push_back(i);
            }
     return indexes;
   }
   
   
   std::vector<int> Statefilter::LocFilter(const Orbitals& orbitals){
     std::vector<int> indexes;
    const std::vector< Eigen::VectorXd >& popE = (_statehist[0].Type()== QMStateType::Singlet)
                   ? orbitals.getFragment_E_localisation_singlet() : orbitals.getFragment_E_localisation_triplet();
    const std::vector< Eigen::VectorXd >& popH = (_statehist[0].Type()== QMStateType::Singlet)
            ? orbitals.getFragment_H_localisation_singlet() : orbitals.getFragment_H_localisation_triplet();
            int fragmentindex=1;
            if (_localiseonA) {
              fragmentindex=0;
            }
            for (unsigned i = 0; i < popE.size(); i++) {
              if (popE[i](fragmentindex) > _loc_threshold && popH[i](fragmentindex) > _loc_threshold) {
                indexes.push_back(i);
              }
            }
     return indexes;
   }
   
   std::vector<int> Statefilter::DeltaQFilter(const Orbitals& orbitals){
    std::vector<int> indexes;
    const std::vector< Eigen::VectorXd >& dQ_frag = (_statehist[0].Type()== QMStateType::Singlet)
                    ? orbitals.getFragmentChargesSingEXC() : orbitals.getFragmentChargesTripEXC();
    for (unsigned i = 0; i < dQ_frag.size(); i++) {
      if (std::abs(dQ_frag[i](0)) > _dQ_threshold) {
        indexes.push_back(i);
      }
    }
     return indexes;
   }

      std::vector<int> Statefilter::OverlapFilter(const Orbitals & orbitals) {
        std::vector<int> indexes;
        if (_statehist.size() <= 1) {
          _laststatecoeff = orbitals.getStateCoefficients(_statehist[0]);
          indexes=std::vector<int>{_statehist[0].Index()};
          return indexes;
        }

        return indexes;
      }
   
   /*
     

          // quasiparticle filter
          if (_has_overlap_filter) {
            if (iter == 0) {

              // One - to - One LIST in 0th iteration
              for (unsigned i = 0; i < orb_iter_input.QPdiagEnergies().size(); i++) {
                state_index.push_back(i);
              }

            } else {

              BasisSet dftbs;
              dftbs.LoadBasisSet(orb_iter_input.getDFTbasis());
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " << orb_iter_input.getDFTbasis() << flush;
              AOBasis dftbasis;
              dftbasis.AOBasisFill(dftbs, orb_iter_input.QMAtoms());
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filled DFT Basis of size " << dftbasis.AOBasisSize() << flush;

              AOOverlap dftoverlap;
              dftoverlap.Fill(dftbasis);
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << dftoverlap.Matrix().rows() << flush;

              // 'LAMBDA' matrix of the present iteration
              Eigen::MatrixXd lambda_N = orb_iter_input.CalculateQParticleAORepresentation();

              // 'LAMBDA' matrix of the previous iteration
              string runFolder_N_1 = jobFolder + "/iter_" + boost::lexical_cast<string>(iter - 1);
              string orbfile_N_1 = runFolder_N_1 + "/system.orb";
              Orbitals orbitals_N_1;
              // load the QM data from serialized orbitals object

              CTP_LOG(ctp::logDEBUG, *_log) << " Loading QM data from " << orbfile_N_1 << flush;
              orbitals_N_1.ReadFromCpt(orbfile_N_1);

              Eigen::MatrixXd lambda_N_1 = orbitals_N_1.CalculateQParticleAORepresentation();
              // calculate QP overlaps
              Eigen::MatrixXd qpoverlaps = lambda_N * dftoverlap.Matrix() * lambda_N_1.transpose();

              // filter for max absolute value (hopefully close to 1)
              for (unsigned j = 0; j < qpoverlaps.cols(); j++) {             
                int maxi = 0;
                double maximumoverlap=qpoverlaps.col(j).maxCoeff(&maxi);
                state_index.push_back(maxi);
                CTP_LOG(ctp::logDEBUG, *_log) << " [" << maxi << " , " << j << "]: " <<maximumoverlap << flush;
              }
            }
          }

          if (_has_osc_filter) {
            // go through list of singlets
            const std::vector<double>oscs = orb_iter_input.Oscillatorstrengths();
            for (unsigned i = 0; i < oscs.size(); i++) {
              if (oscs[i] > _osc_threshold) state_index.push_back(i);
            }
          } else {
            const VectorXfd & energies = (_type == "singlet")
                    ? orb_iter_input.BSESingletEnergies() : orb_iter_input.BSETripletEnergies();
            for (unsigned i = 0; i < energies.size(); i++) {
              state_index.push_back(i);
            }
          }


          // filter according to charge transfer, go through list of excitations in _state_index
          if (_has_dQ_filter) {
            std::vector<int> state_index_copy;
            const std::vector< Eigen::VectorXd >& dQ_frag = (_type == "singlet")
                    ? orb_iter_input.getFragmentChargesSingEXC() : orb_iter_input.getFragmentChargesTripEXC();
            for (unsigned i = 0; i < state_index.size(); i++) {
              if (std::abs(dQ_frag[state_index[i]](0)) > _dQ_threshold) {
                state_index_copy.push_back(state_index[i]);
              }
            }
            state_index = state_index_copy;
          } else if (_has_loc_filter) {
            std::vector<int> state_index_copy;
            const std::vector< Eigen::VectorXd >& popE = (_type == "singlet")
                    ? orb_iter_input.getFragment_E_localisation_singlet() : orb_iter_input.getFragment_E_localisation_triplet();
            const std::vector< Eigen::VectorXd >& popH = (_type == "singlet")
                    ? orb_iter_input.getFragment_H_localisation_singlet() : orb_iter_input.getFragment_H_localisation_triplet();
            if (_localiseonA) {
              for (unsigned i = 0; i < state_index.size(); i++) {
                if (popE[state_index[i]](0) > _loc_threshold && popH[state_index[i]](0) > _loc_threshold) {
                  state_index_copy.push_back(state_index[i]);
                }
              }
            } else {
              for (unsigned i = 0; i < state_index.size(); i++) {
                if (popE[state_index[i]](1) > _loc_threshold && popH[state_index[i]](1) > _loc_threshold) {
                  state_index_copy.push_back(state_index[i]);
                }
              }
            }
            state_index = state_index_copy;
          }


          if (state_index.size() < 1) {
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " WARNING: FILTER yielded no state. Taking lowest excitation" << flush;
            state_index.push_back(0);
          } else {
            if (_type == "quasiparticle") {
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filter yielded QP index: " << state_index[_state - 1 - orb_iter_input.getGWAmin()] << flush;
            } else {
              CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filter yielded state" << _type << ":" << state_index[_state - 1] + 1 << flush;
            }
          }
    */

    }
}
