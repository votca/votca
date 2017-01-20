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

#ifndef __QMMACHINE__H
#define	__QMMACHINE__H

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
// add gwbse header for excited state support
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/gdma.h>


namespace votca { namespace xtp {

    namespace CTP = votca::ctp;
    
// ========================================================================== //
// QM-MM INTERFACE CLASS - CONVERTS BETWEEN QMATOMS <> POLAR OBJECTS          //
// ========================================================================== //
    
class QMMInterface
{
public:
    
    QMMInterface() { _polar_table = CTP::POLAR_TABLE(); };
   ~QMMInterface() {};
    
    // CONVERSION QM -> MM
    CTP::APolarSite *Convert(CTP::QMAtom *atm, int id = -1) {
        double A_to_nm = 0.1;
        CTP::vec pos = A_to_nm*vec(atm->x, atm->y, atm->z);
        double q = atm->charge;
        std::string elem = atm->type;
        double pol = 0.0;
        try {
            pol = _polar_table.at(elem);
        }
        catch(out_of_range) {
            std::cout << std::endl << "QMMInterface - no default polarizability given "
                << "for element type '" << elem << "'. Defaulting to 1A**3" << std::flush;
            pol = 1e-3;
        }

        CTP::APolarSite *new_aps = new CTP::APolarSite(id, elem);
        new_aps->setRank(0);
        new_aps->setPos(pos);
        new_aps->setQ00(q,0); // <- charge state 0 <> 'neutral'
        new_aps->setIsoP(pol);
        
        return new_aps;
    }
    
    CTP::PolarSeg *Convert(std::vector<CTP::QMAtom*> &atms) {        
        CTP::PolarSeg *new_pseg = new CTP::PolarSeg();
        std::vector<CTP::QMAtom*>::iterator it;
        for (it = atms.begin(); it < atms.end(); ++it) {
            CTP::APolarSite *new_site = this->Convert(*it);
            new_pseg->push_back(new_site);
        }
        return new_pseg;
    }
    
    // TODO CONVERSION MM -> QM
    CTP::QMAtom *Convert(CTP::APolarSite*);
    std::vector<CTP::QMAtom*> Convert(CTP::PolarSeg*);
    
private:
    
    // Allocates polarizabilities in A**3 to element types
    std::map<std::string,double> _polar_table;
    
};



// ========================================================================== //
// QM-MM ITERATION CLASS - OBSERVES CONVERGENCE LOOP                          //
// ========================================================================== //
    
class QMMIter
{

public:

    QMMIter(int id) : _id(id), _hasdRdQ(false), _hasQM(false), _hasMM(false)  { ; }
   ~QMMIter() { ; }

   void ConvertPSitesToQMAtoms(std::vector< CTP::PolarSeg* > &, std::vector< CTP::QMAtom* > &);
   void ConvertQMAtomsToPSites(std::vector< CTP::QMAtom* > &, std::vector< CTP::PolarSeg* > &);
   void UpdatePosChrgFromQMAtoms(std::vector< CTP::QMAtom* > &, std::vector< CTP::PolarSeg* > &);  
   void UpdateMPSFromGDMA( std::vector<std::vector<double> > &multipoles,  std::vector< CTP::PolarSeg* > &psegs);
   void GenerateQMAtomsFromPolarSegs(CTP::PolarTop *ptop, Orbitals &orb, bool split_dpl, double dpl_spacing);   

   void setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM);
   void setQMSF(double energy_QM, double energy_SF, double energy_GWBSE);
   void setE_FM(double ef00, double ef01, double ef02, 
                  double ef11, double ef12, double em0,
                  double em1,  double em2, double efm);

   double getRMSdR() { return _dR_RMS; }
   double getRMSdQ() { return _dQ_RMS; }
   double getSUMdQ() { return _dQ_SUM; }
   int getId() { return _id;}

   double getSFEnergy() { assert(_hasQM); return _e_SF; }
   double getFMEnergy() { assert(_hasMM); return _e_fm_; }
   double getQMEnergy() { assert(_hasQM); return _e_QM; }
   double getGWBSEEnergy() { assert(_hasGWBSE); return _e_GWBSE; }
   double getMMEnergy();
   double getQMMMEnergy();


private:

    int    _id;

    bool   _hasdRdQ;
    bool   _hasQM;
    bool   _hasMM;
    bool   _hasGWBSE;

    double _dR_RMS;
    double _dQ_RMS;
    double _dQ_SUM;       

    double _e_QM;
    double _e_SF;
    double _e_GWBSE;

    double _ef_00;
    double _ef_01;
    double _ef_02;
    double _ef_11;
    double _ef_12;
    double _em_0_;
    double _em_1_;
    double _em_2_;
    double _e_fm_;


};
    
    
    
// ========================================================================== //
// QMMACHINE: REGISTER QMPACKAGE TYPE (E.G. GAUSSIAN) AT THE END OF .CC FILE  //
// ========================================================================== //    
    
template< class QMPackage >
class QMMachine 
{
    
public:

    QMMachine(CTP::XJob *job, CTP::XInductor *xind, QMPackage *qmpack,
              CTP::Property *opt, string sfx, int nst, bool mav);
   ~QMMachine();
    
    void Evaluate(CTP::XJob *job);
    //void WriteQMPackInputFile(string inputFile, QMPackage *qmpack, XJob *job);
    
    bool Iterate(string jobFolder, int iterCnt);    
    QMMIter *CreateNewIter();
    bool hasConverged();
    bool AssertConvergence() { return _isConverged; }
    
    void setLog(CTP::Logger *log) { _log = log; }
    
private:    
    
    CTP::XJob *_job;
    CTP::XInductor *_xind;
    QMPackage *_qmpack;
    CTP::Logger *_log;
    int _subthreads;
    
    std::vector<QMMIter*> _iters;
    bool _isConverged;
    int _maxIter;

    // GDMA object
    // GDMA _gdma;
    Property _gdma_options;
    bool _do_gdma;
    
    
    
    
    
    // GWBSE object
    // GWBSE _gwbse;
    Property _gwbse_options;
    int      _state;
    string   _type;
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
    
    bool _do_gwbse; // needs to be set by options!!!

};


}}

#endif
