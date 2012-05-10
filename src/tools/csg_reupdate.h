/*
 * File:   csg_reupdate.h
 * Author: mashaya1
 *
 * Created on October 13, 2011, 11:09 PM
 */

#ifndef CSG_REUPDATE_H
#define	CSG_REUPDATE_H
#include <boost/program_options.hpp>
#include <votca/csg/csgapplication.h>
#include <votca/tools/table.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/property.h>
#include <votca/tools/histogramnew.h>
#include "potentialfunctions/potentialfunction.h"
#include "potentialfunctions/potentialfunctioncbspl.h"
#include "potentialfunctions/potentialfunctionlj126.h"
#include "potentialfunctions/potentialfunctionljg.h"
#include <votca/csg/topologyreader.h>

using namespace votca::csg;
using namespace votca::tools;
using namespace std;

struct PotentialInfo {

        PotentialInfo(int index,
		      bool bonded_,
		      int vec_pos_, 
		      string& param_in_ext_, string& rdf_ext_,
		      Property *options);

        int potentialIndex;
        bool bonded;
        PotentialFunction *ucg;
        int vec_pos;
        pair<int, int> beadTypes;

        string potentialName;
        string potentialFunction;
        string type1, type2;

        Table aardf;
        double rdf_norm;

        double rmin,rcut;

        Property *_options;
};

class CsgREupdate
: public CsgApplication
{
public:
    string ProgramName() { return "csg_reupdate"; }
    void HelpText(ostream &out) {
        out << "computes relative entropy update.";
            }

    bool DoTrajectory() { return true;}
  
    bool DoMapping(){ return false; }

    bool DoThreaded() { return true; }
    bool SynchronizeThreads() { return false; }

    bool NeedsTopology() { return false; }

    void Initialize();
    bool EvaluateOptions();
    void BeginEvaluate(Topology *top, Topology *top_atom = 0);
    void LoadOptions(const string &file);
    
    void Run();
    
    void EndEvaluate();
    CsgApplication::Worker *ForkWorker(void);
    void MergeWorker(Worker *worker);
    
private:

protected:
   
    Property _options;
    list<Property *> _nonbonded;

    typedef vector<PotentialInfo *> PotentialContainer;
    PotentialContainer _potentials;

    int _nlamda;
    ub::vector<double> _lamda;
    ub::vector<double> _dlamda;
    // _HS is a symmetric matrix
    ub::symmetric_matrix<double, ub::upper> _HS;
    ub::vector<double> _DS;
    ub::vector<double> _dUFrame;
    
    double _UavgAA;
    double _UavgCG;
    double _beta;
    double _relax;
    int _nframes;
   
    bool _gentable;
    bool _dosteep;
    
    // file extension for the inputs/outputs
    string _param_in_ext, _param_out_ext;
    string _pot_out_ext;
    string _rdf_ext;

    void WriteOutFiles();
    void EvalBonded(Topology *conf, PotentialInfo *potinfo);
    void EvalNonbonded(Topology *conf, PotentialInfo *potinfo);

    // Compute Avg U, dU, and d2U values in reference AA ensemble
    void AAavgBonded(PotentialInfo *potinfo);
    void AAavgNonbonded(PotentialInfo *potinfo);

    // Formulates _HS dlamda = - _DS system of Lin Eq.
    void REFormulateLinEq();

    // Solve _HS dlamda = - _DS and update _lamda
    void REUpdateLamda();
   
};


class CsgREupdateWorker
: public CsgApplication::Worker
{
public:

    ~CsgREupdateWorker(){};

    Property _options;
    list<Property *> _nonbonded;
    
    typedef vector<PotentialInfo *> PotentialContainer;
    PotentialContainer _potentials;

    int _nlamda;
    ub::vector<double> _lamda;
    ub::symmetric_matrix<double, ub::upper> _HS;
    ub::vector<double> _DS;
    ub::vector<double> _dUFrame;

    double _UavgAA;
    double _UavgCG;
    double _beta;
    int _nframes;

    void EvalConfiguration(Topology *conf, Topology *conf_atom);
    void EvalBonded(Topology *conf, PotentialInfo *potinfo);
    void EvalNonbonded(Topology *conf, PotentialInfo *potinfo);

};

#endif	/* CSG_REUPDATE_H */
