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

#ifndef VOTCA_XTP_QMAPEMACHINE_H
#define	VOTCA_XTP_QMAPEMACHINE_H



#include <votca/xtp/dftengine.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/espfit.h>

#include <votca/xtp/ewaldnd.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/xinductor.h>

#include <votca/xtp/qminterface.h>
#include <votca/xtp/qmiter.h>
#include <votca/xtp/statefilter.h>

namespace votca { namespace xtp {



class QMAPEMachine
{
    
public:

    QMAPEMachine(xtp::XJob *job, xtp::Ewald3DnD *cape, 
              Property *opt, std::string sfx);
   ~QMAPEMachine();
    
    void Evaluate(xtp::XJob *job);
 
    bool AssertConvergence() { return _isConverged; }
    
    void setLog(xtp::Logger *log) { _log = log; }
    
private:
    QMMIter *CreateNewIter();
    bool EvaluateGWBSE(Orbitals &orb, std::string runFolder);
    bool hasConverged();
    bool Iterate(std::string jobFolder, int iterCnt);
    void SetupPolarSiteGrids( const std::vector< const vec *>& gridpoints,const std::vector< QMAtom* >& atoms);
    std::vector<double> ExtractNucGrid_fromPolarsites();
    std::vector<double> ExtractElGrid_fromPolarsites();
    
    Statefilter _filter;
    QMInterface qminterface;
    xtp::Logger *_log;

    bool _run_ape;
    bool _run_dft;
    bool _run_gwbse;

    xtp::XJob *_job;
    xtp::Ewald3DnD *_cape;

    std::vector<QMMIter*> _iters;
    int _maxIter;
    bool _isConverged;
    
    unsigned NumberofAtoms;
    string _externalgridaccuracy;
    
    Property _gwbse_options;
    Property _dft_options;
    QMState  _initialstate;
    
    double _crit_dR;
    double _crit_dQ;
    double _crit_dE_QM;
    double _crit_dE_MM;   
    
    std::vector< xtp::PolarSeg* > target_bg;     
    std::vector< xtp::PolarSeg* > target_fg;     

};


}}

#endif // VOTCA_XTP_QMAPEMACHINE_H

