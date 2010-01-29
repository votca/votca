#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>
#include <stdexcept>
#include "easyjobserver.h"
#include "qmtopology.h"

using namespace std;

void help_text(void)
{
    votca::csg::HelpTextHeader("easyJ");
    cout << "Calculate transfer integrals\n\n";
}

 /// Namespace to read in program options
namespace po=boost::program_options;

void check_option(po::options_description &desc, po::variables_map &vm, const string &option)
{
    if(!vm.count(option)) {
        cout << "easyJ \n\n";
        cout << desc << endl << "parameter " << option << " is not specified\n";
        exit(1);
    }
}

int main(int argc, char** argv)
{    
    int write_every=0;
    // we have one observer
    EasyJObserver observer;
    Property options;
    // The CGEngine does the work
    CGEngine cg_engine;
    QMTopology qmtopol;

    try {
                // let cg_engine add some program options
        cg_engine.Initialize();
   
        cg_engine.AddProgramOptions()
            ("listcharges,l", po::value<string>(), "  Crg unit definitions");
        cg_engine.AddProgramOptions()
            ("cutoff,c", po::value<double>()-> default_value(1.0), "  CutOff for nearest neighbours");
        cg_engine.AddProgramOptions()
            ("nnnames,n", po::value<string>()-> default_value("*"), "  List of strings that the concatenation of the two molnames must match to be printed");
        /// Parameters required to calculate rates and to run KMC
        cg_engine.AddProgramOptions()
            ("options,o", po::value<string>(), "  KMC and MD2QM options");
        cg_engine.ParseCommandLine(argc, argv);

        po::variables_map &vm
            = cg_engine.OptionsMap();
    
        // Help menu
        if (vm.count("help")) {
            help_text();
            cout << cg_engine.OptionsDesc() << endl;
            return 0;
        }

        check_option(desc, cg_engine.OptionsMap(), "options");

        load_property_from_xml(options, vm["options"].as<string>());
        observer.Initialize(qmtopol, options);
        // add our observer that it gets called to analyze frames
        cg_engine.AddObserver((CGObserver*)&observer);

        qmtopol.LoadListCharges(vm["listcharges"].as<string>());
        observer.setCutoff(vm["cutoff"].as<double>());
        observer.setNNnames(vm["nnnames"].as<string>());
        // try to run the cg process, go through the frames, etc...
        cg_engine.Run();

    }

    /// error handling
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
        return -1;
    }
    return 0;
}


