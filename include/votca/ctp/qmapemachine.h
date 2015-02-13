#ifndef __QMAPEMACHINE__H
#define	__QMAPEMACHINE__H


#include <votca/ctp/votca_ctp_config.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/gwbse.h>
#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/espfit.h>
#include <votca/ctp/ewaldnd.h>



namespace votca { namespace ctp {

    
// ========================================================================== //
// QM-MM INTERFACE CLASS - CONVERTS BETWEEN QMATOMS <> POLAR OBJECTS          //
// ========================================================================== //
    
class QMAPEInterface
{
public:
    
    QMAPEInterface() { _polar_table = POLAR_TABLE(); };
   ~QMAPEInterface() {};
    
    // CONVERSION QM -> MM
    APolarSite *Convert(QMAtom *atm, int id = -1) {
        double A_to_nm = 0.1;
        vec pos = A_to_nm*vec(atm->x, atm->y, atm->z);
        double q = atm->charge;
        string elem = atm->type;
        double pol = 0.0;
        try {
            pol = _polar_table.at(elem);
        }
        catch(out_of_range) {
            cout << endl << "QMAPEInterface - no default polarizability given "
                << "for element type '" << elem << "'. Defaulting to 1A**3" << flush;
            pol = 1e-3;
        }

        APolarSite *new_aps = new APolarSite(id, elem);
        new_aps->setRank(0);
        new_aps->setPos(pos);
        new_aps->setQ00(q,0); // <- charge state 0 <> 'neutral'
        new_aps->setIsoP(pol);
        
        return new_aps;
    }
    
    PolarSeg *Convert(vector<QMAtom*> &atms) {        
        PolarSeg *new_pseg = new PolarSeg();
        vector<QMAtom*>::iterator it;
        for (it = atms.begin(); it < atms.end(); ++it) {
            APolarSite *new_site = this->Convert(*it);
            new_pseg->push_back(new_site);
        }
        return new_pseg;
    }
    
    // TODO CONVERSION MM -> QM
    QMAtom *Convert(APolarSite*);
    vector<QMAtom*> Convert(PolarSeg*);
    
private:
    
    // Allocates polarizabilities in A**3 to element types
    map<string,double> _polar_table;
    
};



// ========================================================================== //
// QM-MM ITERATION CLASS - OBSERVES CONVERGENCE LOOP                          //
// ========================================================================== //
    
class QMAPEIter
{

public:

    QMAPEIter(int id) : _hasdRdQ(false), _hasQM(false), _hasMM(false), _id(id) { ; }
   ~QMAPEIter() { ; }

   void ConvertPSitesToQMAtoms(vector< PolarSeg* > &, vector< QMAtom* > &);
   void ConvertQMAtomsToPSites(vector< QMAtom* > &, vector< PolarSeg* > &);
   void UpdatePosChrgFromQMAtoms(vector< QMAtom* > &, vector< PolarSeg* > &);   
   void GenerateQMAtomsFromPolarSegs(vector<PolarSeg*> &qm, vector<PolarSeg*> &mm, Orbitals &orb);

   void setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM);
   void setQMSF(double energy_QM, double energy_SF, double energy_GWBSE);
   void setE_FM(double ef00, double ef01, double ef02, 
                  double ef11, double ef12, double em0,
                  double em1,  double em2, double efm);

   double getRMSdR() { return _dR_RMS; }
   double getRMSdQ() { return _dQ_RMS; }
   double getSUMdQ() { return _dQ_SUM; }

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
class QMAPEMachine
{
    
public:

	QMAPEMachine(XJob *job, Ewald3DnD *cape, QMPackage *qmpack,
              Property *opt, string sfx, int nst);
   ~QMAPEMachine();
    
    void Evaluate(XJob *job);
    //void WriteQMPackInputFile(string inputFile, QMPackage *qmpack, XJob *job);
    
    bool Iterate(string jobFolder, int iterCnt);
    bool EvaluateGWBSE(Orbitals &orb, string runFolder);
    QMAPEIter *CreateNewIter();
    bool hasConverged();
    bool AssertConvergence() { return _isConverged; }
    
    void setLog(Logger *log) { _log = log; }
    
private:    
    

    Logger *_log;
    int _subthreads;

    bool _run_ape;
    bool _run_dft;
    bool _run_gwbse;

    XJob *_job;
    XInductor *_xind;
    QMPackage *_qmpack;
    Ewald3DnD *_cape;

    vector<QMAPEIter*> _iters;
    int _maxIter;
    bool _isConverged;

    // GWBSE object
    GWBSE _gwbse;
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

};


}}

#endif
