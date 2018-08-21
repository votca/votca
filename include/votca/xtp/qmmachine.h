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

#ifndef __QMMACHINE__H
#define	__QMMACHINE__H


#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>

// add gwbse header for excited state support
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/gdma.h>
#include <votca/xtp/qminterface.h>
#include <votca/xtp/qmiter.h>

namespace votca { namespace xtp {



class QMMachine{

public:

    QMMachine(ctp::XJob *job, ctp::XInductor *xind, QMPackage *qmpack,
              Property *opt, string sfx, int nst, bool mav);
   ~QMMachine();

    int Evaluate(ctp::XJob *job);

    void setLog(ctp::Logger *log) { _log = log; }

private:
    bool Iterate(string jobFolder, int iterCnt);
    bool RunDFT(string& runFolder, std::vector<std::shared_ptr<ctp::PolarSeg> >& MultipolesBackground);
    void RunGWBSE(string& runFolder);
    void RunGDMA(QMMIter* thisIter, string& runFolder);
    void Density2Charges( std::vector<int> state_index ={});
    
    QMMIter *CreateNewIter();
    bool hasConverged();
    ctp::XJob *_job;
    ctp::XInductor *_xind;
    QMPackage *_qmpack;
    ctp::Logger *_log;
    int _subthreads;

    std::vector<QMMIter*> _iters;
    bool _isConverged;
    int _maxIter;

    Property _gdma_options;
    bool _do_gdma;


    Property _gwbse_options;
    int      _state;
    string   _type;
    bool     _has_osc_filter=false;
    bool     _has_overlap_filter=false;
    double   _osc_threshold;
    bool     _has_dQ_filter=false;
    bool     _has_loc_filter=false;
    double   _dQ_threshold;
    double   _loc_threshold;
    bool     _localiseonA=false;
    double _crit_dR;
    double _crit_dQ;
    double _crit_dE_QM;
    double _crit_dE_MM;

    bool _do_gwbse; // needs to be set by options!!!
    bool _do_archive;
    bool _static_qmmm;
    Orbitals orb_iter_input;
    
    double _alpha;
    Eigen::MatrixXd _DMAT_old;

    

};


}}

#endif
