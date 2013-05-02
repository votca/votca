#ifndef __QMMACHINE__H
#define	__QMMACHINE__H


#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/gaussian.h>
#include <votca/ctp/orbitals.h>


namespace votca { namespace ctp {


class QMMachine 
{

public:

    QMMachine(XJob *job, XInductor *xind, Gaussian *qmpack,
              Property *opt, string sfx, int nst, bool mav);
   ~QMMachine();
   
    // ================================================= //
    // QM-MM ITERATION CLASS - OBSERVES CONVERGENCE LOOP //
    // ================================================= //
   
    class QMMIter 
    {
        
    public:
        
        QMMIter(int id) : _hasdRdQ(false), _hasQM(false), _hasMM(false), _id(id) { ; }
       ~QMMIter() { ; }
        
       void ConvertPSitesToQMAtoms(vector< PolarSeg* > &, vector< QMAtom* > &);
       void ConvertQMAtomsToPSites(vector< QMAtom* > &, vector< PolarSeg* > &);
       void UpdatePosChrgFromQMAtoms(vector< QMAtom* > &, vector< PolarSeg* > &);       
       
       void setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM);
       void setQMSF(double energy_QM, double energy_SF);
       void setE_FM(double ef00, double ef01, double ef02, 
                      double ef11, double ef12, double em0,
                      double em1,  double em2, double efm);
       
       double getRMSdR() { return _dR_RMS; }
       double getRMSdQ() { return _dQ_RMS; }
       double getSUMdQ() { return _dQ_SUM; }
       
       double getSFEnergy() { assert(_hasQM); return _e_SF; }
       double getFMEnergy() { assert(_hasMM); return _e_fm_; }
       double getQMEnergy() { assert(_hasQM); return _e_QM; }
       double getMMEnergy();
       double getQMMMEnergy();
        
        
    private:
        
        int    _id;
        
        bool   _hasdRdQ;
        bool   _hasQM;
        bool   _hasMM;
        
        double _dR_RMS;
        double _dQ_RMS;
        double _dQ_SUM;       
        
        double _e_QM;
        double _e_SF;        
        
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

    // ================================================= //
    // EVALUATE - ITERATE - TEST CONVERGENCE             //
    // ================================================= //
    
    void Evaluate(XJob *job);
    void WriteQMPackInputFile(string inputFile, Gaussian *qmpack, XJob *job);
    bool Iterate(string jobFolder, int iterCnt);
    QMMIter *CreateNewIter();
    bool hasConverged();
    bool AssertConvergence() { return _isConverged; }
    
private:

    XJob *_job;
    XInductor *_xind;
    Gaussian *_qmpack;
    int _subthreads;
    
    vector<QMMIter*> _iters;
    bool _isConverged;
    int _maxIter;
    
    double _crit_dR;
    double _crit_dQ;
    double _crit_dE_QM;
    double _crit_dE_MM;
    
    bool _convg_dR;
    bool _convg_dQ;
    bool _convg_dE_QM;
    bool _convg_dE_MM;

};

    
}}

#endif