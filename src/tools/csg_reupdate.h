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
#include <boost/numeric/ublas/matrix.hpp>
#include "potfunctions.h"
#include <votca/tools/property.h>
#include "typedefpotfun.h"

using namespace votca::csg;
using namespace votca::tools;
using namespace std;


class CsgREupdate
: public CsgApplication
{
public:
    string ProgramName() { return "csg_reupdate"; }
    void HelpText(ostream &out) {
        out << "computes relative entropy update.";
            }

    bool DoTrajectory() {return true;}
    // we do not want to perform mapping
    bool DoMapping() { return false;}

    void Initialize();
    bool EvaluateOptions();
    void BeginEvaluate(Topology *top, Topology *top_atom);
    //void Run();
    // write out results in EndEvaluate
    void EndEvaluate();
    // do RE update calculation in this function
    void EvalConfiguration(Topology *conf, Topology *conf_atom);
    void LoadOptions(const string &file);
    
private:

protected:
    // define what functional form to use for CG potential
    // defined in typedefpotfun.h
    //typedef FunctionLJ126 PotFunction;
    
    struct PotentialInfo {
        PotentialInfo(int index, bool bonded_, int vec_pos_, Property *options);
        int potentialIndex;
        bool bonded;
        PotFunction ucg;
        int vec_pos;
        pair<int, int> beadTypes;
        
        string potentialName;
        string type1, type2;

        Table aahist;

        ub::vector<double> pottblgrid;


        Property *_options;
    };
    Property _options;
    list<Property *> _nonbonded;

    typedef vector<PotentialInfo *> PotentialContainer;
    PotentialContainer _potentials;

    int _nlamda;
    ub::vector<double> _lamda;
    ub::matrix<double> _HS;
    ub::vector<double> _DS;

    double _beta;
    double _relax;
    int _nframes;


    void WriteOutFiles();
    void REFormulateLinEq();
    void REUpdateLamda();
    void EvalBonded(Topology *conf, PotentialInfo *potinfo);
    void EvalNonbonded(Topology *conf, PotentialInfo *potinfo);
    void AAavgBonded(PotentialInfo *potinfo);
    void AAavgNonbonded(PotentialInfo *potinfo);
    
};
#endif	/* CSG_REUPDATE_H */
