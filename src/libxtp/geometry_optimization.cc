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

#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/forces.h>
#include <votca/xtp/bfgs-trm.h>

#include "votca/xtp/statefilter.h"

namespace votca {
    namespace xtp {

        void GeometryOptimization::Initialize(tools::Property &options) {
           
            std::string statestring= options.ifExistsReturnElseThrowRuntimeError<std::string>(".state");
            _opt_state.FromString(statestring);
            if(!_opt_state.Type().isExciton()){
              throw std::runtime_error("At the moment only excitonic states can be optimized");
            }
            std::vector<std::string> choices = {"BFGS-TRM"};
            _optimizer = options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(".optimizer.method", choices);
            _optimizer_options = options.get(".optimizer");

            if(options.exists(".forces")){
            _force_options = options.get(".forces");
            }else{
                throw std::runtime_error("No forces options provided");
            }
            if(options.exists(".filter")){
                _filter_options=options.get(".filter");
            }else{
                throw std::runtime_error("No filter options set");
            }

            return;
        }

        void GeometryOptimization::Evaluate() {
            CTP_LOG(ctp::logINFO, *_pLog) << "Requested geometry optimization of excited state " <<  _opt_state.ToString() << std::flush;
            
            Statefilter filter;
            filter.Initialize(_filter_options);
            filter.setInitialState(_opt_state);
            filter.setLogger(_pLog);
            filter.PrintInfo();

            // get a force object
            Forces force_engine(_gwbse_engine, _qmpackage,_orbitals);
            force_engine.Initialize(_force_options);
            force_engine.setLog(_pLog);

            // get the optimizer
            if (_optimizer == "BFGS-TRM") {
                BFGSTRM bfgstrm(_gwbse_engine, _qmpackage, _orbitals, force_engine);
                bfgstrm.Initialize(_optimizer_options);
                bfgstrm.setLog(_pLog);
                bfgstrm.Optimize();
            }
            return;
        }

    }
}
