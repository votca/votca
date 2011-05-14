/* 
 * File:   cont_kmc.h
 * Author: vehoff
 *
 * Created on April 8, 2010, 12:04 PM
 */

#ifndef _CONT_KMC_H
#define	_CONT_KMC_H

#include "qmcalculator.h"
#include <votca/kmcold/vertex.h>
#include <votca/kmcold/hoppers.h>
#include <votca/kmcold/kmc.h>
#include <votca/kmcold/graph.h>
#include <time.h>

class ContKmc : public QMCalculator{
public:
    ContKmc() {};
    ~ContKmc() {};

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);

private:
    /// electric field
    vec _E;
    /// total simulation time
    double _total_time;
    /// constant
    double _alpha;
    /// time step
    double _dt;
    /// number of charges simultaneously present in a simultion
    int _ncrg;
    /// number of KMC runs to be performed
    int _nruns;
    ///  output streams for velocity averaging & diffusion
    ofstream _out_cont;
    ofstream _out_diff;
    unsigned int _seed;
    /// creation of KMC graph
    void make_kmc_graph(QMTopology *top, graph *a, QMNBList &nblist);
};

inline void ContKmc::Initialize(QMTopology *top, Property *options){
    /// Read in the electric field - should be the same as used for the rates
    _E = options->get("options.calc_rates.e_field").as<vec>();

    /// Read simulation parameters from input file (main.xml)
    _total_time = options->get("options.kmc_cont.total_time").as<double>();
    _dt = options->get("options.kmc_cont.dt").as<double>();
    _alpha = options->get("options.kmc_cont.alpha").as<double>();
    _ncrg = options->get("options.kmc_cont.ncrg").as<int>();
    _nruns = options->get("options.kmc_cont.nruns").as<int>();

    _seed = time(NULL)%0xFFFFFFFF;
    if(options->exists("options.kmc_cont.seed"))
        _seed = options->get("options.kmc_cont.seed").as<unsigned int>();

    /// Initialize output files for continuous KMC and diffusion
    _out_cont.open("kmc_cont.res");
    _out_diff.open("kmc_diff.res");
    if(_out_cont!=0){
        _out_cont << "# seed for random number generator: " << _seed << endl;
        _out_cont << "# number of run, time [sec], average velocity [m/s], mobility [cm^2/sec], electric field [V/m], end to end vector [nm]" << endl;
    }
    if(_out_diff!=0){
        _out_diff << "# time [sec], displacement in x, y and z direction [m]" << endl;
    }
    _out_cont.close();
    _out_diff.close();

    /// Initialize the random number generator
    // \todo this routine is shitty, replace it!
    // according to marsaglia paper, choosing 4 initial values out of range
    // demanded by the generator will
    // give satisfactory results, however doesn't guarantee same execution on
    // different machins. This also applies for rand() which might be specific
    // for different architectures.
    srand(_seed);
    Random::init(rand(), rand(), rand(), rand());
}

inline bool ContKmc::EvaluateFrame(QMTopology *top){
    /// creating graph and initializing generators and hoppers
    QMNBList &nblist = top->nblist();
    graph kmc_grid;
    make_kmc_graph(top,&kmc_grid,nblist);
    kmc_grid.setGeneratorsOnly();
    hoppers charges(&kmc_grid);
    charges.setRecordOcc(true);

    /// preparing KMC Algorithm
    KMCAlg kmc_alg(_total_time, _dt, _alpha, _ncrg, &charges);
    cout << "Starting continuos KMC." << endl;

    /// A single KMC run
    kmc_alg.kmcPBC(_nruns, _out_cont, _out_diff);

    // copy back occupation probabilities
    vector<vertex*>::iterator it = kmc_grid.getFirstVertex();

    for (; it != kmc_grid.getLastVertex(); ++it) {
        QMCrgUnit *crg = top->GetCrgUnit((*it)->getCrgUnitId());
        if(crg == NULL)
            throw std::runtime_error("did not find crgunit id given in kmc graph");
        crg->setOccupationProbability(charges.getOccProbability(*it));
    }


    cout << "Finished continuous KMC." << endl;
}

inline void ContKmc::make_kmc_graph(QMTopology *top, graph *a, QMNBList &nblist) {
    cout << "[make_kmc_graph]: Building KMC Graph...";
    /// assign constants
    a->SetField(_E);
    /// set vertices equal to centers of mass
    vector < QMCrgUnit *> &listCharges = top->CrgUnits();
    vector < QMCrgUnit *>::iterator it;
    for (it = listCharges.begin(); it != listCharges.end(); ++it) {
        vertex *v = a->AddVertex((*it)->GetCom(), _E); /// TO DO: remove necessity for E-field at this point
        v->setCrgUnitId((*it)->getId());
    }
    /// set edges, two edges 1->2 and 2->1 are created per neighboring pair
    for (QMNBList::iterator iter = nblist.begin(); iter != nblist.end(); ++iter) {
        a->AddEdge((*iter)->first->getId(), (*iter)->second->getId(), (*iter)->rate12(), (*iter)->r());
        a->AddEdge((*iter)->second->getId(), (*iter)->first->getId(), (*iter)->rate21(), -(*iter)->r());
    }
    cout << " Done." << endl;
}

#endif	/* _CONT_KMC_H */

