#include "kmc_cont_app.h"

KmcCont::KmcCont()
{}

KmcCont::~KmcCont()
{}

void KmcCont::HelpText(){
    cout << "Continuous Kinetic Monte Carlo \n\n";
    cout << _op_desc_specific << endl;
}

void KmcCont::Initialize(){
    /// Read in the electric field - should be the same as used for the rates
    _E = _options.get("options.calc_rates.e_field").as<vec>();

    /// Read simulation parameters from input file (main.xml)
    _total_time = _options.get("options.kmc_cont.total_time").as<double>();
    _dt = _options.get("options.kmc_cont.dt").as<double>();
    _alpha = _options.get("options.kmc_cont.alpha").as<double>();
    _ncrg = _options.get("options.kmc_cont.ncrg").as<int>();
    _nruns = _options.get("options.kmc_cont.nruns").as<int>();

    /// Initialize output files for continuous KMC and diffusion
    _out_cont.open("kmc_cont.res");
    _out_diff.open("kmc_diff.res");
    if(_out_cont!=0){
        _out_cont << "# number of run, time [sec], average velocity [m/s], mobility [cm^2/sec], electric field [V/m], end to end vector [nm]" << endl;
    }
    if(_out_diff!=0){
        _out_diff << "# time [sec], displacement in x, y and z direction [m]" << endl;
    }
    _out_cont.close();
    _out_diff.close();

    /// Initialize the random number generator
    Random::init(14, 122, 472, 1912);
}

void KmcCont::EvaluateFrame(){

    /// creating graph and initializing generators and hoppers
    QMNBList &nblist = _qmtop.nblist();
    graph kmc_grid;
    make_kmc_graph(&kmc_grid,nblist);
    kmc_grid.setGeneratorsOnly();
    hoppers charges(&kmc_grid);
    
    /// preparing KMC Algorithm
    KMCAlg kmc_alg(_total_time, _dt, _alpha, _ncrg, &charges);
    cout << "Starting continuos KMC." << endl;

    /// A single KMC run
    kmc_alg.kmcPBC(_nruns, _out_cont, _out_diff);
    cout << "Finished continuous KMC." << endl;
}

void KmcCont::make_kmc_graph(graph *a, QMNBList &nblist) {
    cout << "[make_kmc_graph]: Building KMC Graph...";
    /// assign constants
    a->SetField(_E);
    /// set vertices equal to centers of mass
    list < CrgUnit *> listCharges = _qmtop.crglist();
    list < CrgUnit *>::iterator it;
    for (it = listCharges.begin(); it != listCharges.end(); ++it) {
        a->AddVertex((*it)->GetCom(), _E); /// TO DO: remove necessity for E-field at this point
    }
    /// set edges, two edges 1->2 and 2->1 are created per neighboring pair
    for (QMNBList::iterator iter = nblist.begin(); iter != nblist.end(); ++iter) {
        a->AddEdge((*iter)->first->getId(), (*iter)->second->getId(), (*iter)->rate12(), (*iter)->r());
        a->AddEdge((*iter)->second->getId(), (*iter)->first->getId(), (*iter)->rate21(), -(*iter)->r());
    }
    cout << " Done." << endl;
}