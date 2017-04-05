/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef __QMAPEMACHINE__H
#define	__QMAPEMACHINE__H




#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/espfit.h>

#include <votca/ctp/ewaldnd.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>

#include <votca/xtp/qminterface.h>
#include <votca/xtp/qmiter.h>

namespace votca { namespace xtp {


    
// ========================================================================== //
// QMMACHINE: REGISTER QMPACKAGE TYPE (E.G. GAUSSIAN) AT THE END OF .CC FILE  //
// ========================================================================== //    
    
template< class QMPackage >
class QMAPEMachine
{
    
public:

	QMAPEMachine(ctp::XJob *job, ctp::Ewald3DnD *cape, QMPackage *qmpack,
              Property *opt, std::string sfx, int nst);
   ~QMAPEMachine();
    
    void Evaluate(ctp::XJob *job);
    bool Iterate(std::string jobFolder, int iterCnt);
    bool EvaluateGWBSE(Orbitals &orb, std::string runFolder);
    QMMIter *CreateNewIter();
    bool hasConverged();
    bool AssertConvergence() { return _isConverged; }
    
    void setLog(ctp::Logger *log) { _log = log; }
    
private:    
    

    ctp::Logger *_log;
    int _subthreads;

    bool _run_ape;
    bool _run_dft;
    bool _run_gwbse;

    ctp::XJob *_job;
    ctp::XInductor *_xind;
    QMPackage *_qmpack;
    ctp::Ewald3DnD *_cape;
    
    

    std::vector<QMMIter*> _iters;
    int _maxIter;
    bool _isConverged;

    // GWBSE object
    
    Property _gwbse_options;
    int      _state;
    std::string   _type;
    bool     _has_osc_filter;
    double   _osc_threshold;
    bool     _has_dQ_filter;
    double   _dQ_threshold;   
    
    double _crit_dR;
    double _crit_dQ;
    double _crit_dE_QM;
    double _crit_dE_MM;
    
    bool _convg_dR;
    bool _convg_dQ;
    bool _convg_dE_QM;
    bool _convg_dE_MM;
    
    bool _split_dpl;
    double _dpl_spacing;
    
    bool   _exportgridtofile;

};


}}

#endif
