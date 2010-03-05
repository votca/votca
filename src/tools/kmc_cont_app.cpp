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

bool KmcCont::EvaluateFrame(){

    /// creating graph and initializing generators and hoppers
    QMNBList &nblist = _qmtop.nblist();
    graph kmc_grid;
    make_kmc_graph(&kmc_grid,nblist);
    kmc_grid.setGeneratorsOnly();
    hoppers charges(&kmc_grid);
    
    /// preparing KMC Algorithm
    KMCAlg kmc_alg(1E-8, 1E-15, 1.2, 1, &charges);
    cout << "Starting continuos KMC." << endl;

    /// A single KMC run
    kmc_alg.kmcPBC(_nruns, _out_cont, _out_diff);
    cout << "Finished continuous KMC." << endl;
}