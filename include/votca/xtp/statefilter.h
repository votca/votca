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

#ifndef _VOTCA_XTP_STATEFILTER_H
#define _VOTCA_XTP_STATEFILTER_H

#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>
#include <votca/xtp/qmstate.h>


namespace votca {
namespace xtp {
/**
 *  \brief  Filters from a spectrum of states the state, which fullfills certain criteria
 *
 *
 */

class Statefilter {

public:
    void Initialize(tools::Property& options);
    void setLogger(ctp::Logger* log){_log=log;}
    void setInitialState(const QMState& state ){_statehist.push_back(state);}
    void PrintInfo()const;
    void Filter(Orbitals& orbitals);
    const QMState& getState(){return _state;}// zero indexed;
    
private:
 
    std::vector<int> OscFilter(const Orbitals& orbitals);
    std::vector<int> LocFilter(const Orbitals& orbitals);
    std::vector<int> DeltaQFilter(const Orbitals& orbitals);
    std::vector<int> OverlapFilter(Orbitals& orbitals);
    
    Eigen::VectorXd CalculateOverlap(Orbitals & orbitals);
    
    void UpdateLastCoeff(Orbitals& orbitals);
    Eigen::MatrixXd CalcOrthoCoeffs(Orbitals& orbitals);

    std::vector<int> CollapseResults(std::vector< std::vector<int> >& results)const;
    std::vector<int> ComparePairofVectors( std::vector<int>& vec1, std::vector<int>& vec2)const;

QMState _initial_state;
QMState _state;    
ctp::Logger *_log;
 
std::vector<QMState> _statehist;

Eigen::VectorXd _laststatecoeff;
Eigen::MatrixXd _S_onehalf;

bool _use_oscfilter=false;
double _oscthreshold=0.0;

bool _use_overlapfilter=false;
double _overlapthreshold;

bool _use_localisationfilter=false;
bool _localiseonA;
double _loc_threshold;

bool _use_dQfilter=false;
double _dQ_threshold=true;


};
}
}

#endif /* _VOTCA_XTP_STATEFILTER_H */
