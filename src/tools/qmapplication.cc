#include "qmapplication.h"

void QMApplication::ReadInput(int argc, char **argv)
{
    namespace po = boost::program_options;

    /// define standard program options
    _op_desc.add_options()
    ("help", "  produce this help message")
    ("o, opt", boost::program_options::value<string>()->default_value("main.xml"), "  main program options")
    ;

    /// define options related to the trajectory
    po::options_description trj("Trajectory options");
    trj.add_options()
    ("begin", boost::program_options::value<double>(), "  skip frames before this time")
    ("first-frame", boost::program_options::value<int>()->default_value(0), "  start with this frame")
    ("nframes", boost::program_options::value<int>()->default_value(1), "  process so many frames")
    ;
    _op_desc.add(trj);

    /// add specific options defined via Initialize of the child class
    _op_desc.add(_op_desc_specific);

    /// parse the command line
    try {
        po::store(po::parse_command_line(argc, argv, _op_desc), _op_vm);
        po::notify(_op_vm);
    }
    catch(boost::program_options::error err) {
        throw runtime_error(string("error parsing command line: ") + err.what());
    }
}

void QMApplication::Run(int argc, char **argv)
{

    HelpText(); /// print helpful information about the child
    Initialize(); /// initialize program-specific parameters
    ReadInput(argc, argv); /// initialize general parameters & read input file
    
    bool has_begin=false; /// was a starting time specified?
    double begin;   /// strating time
    int first_frame; /// starting frame
    int nframes = 1; /// number of frames to be processed


    if (!_op_vm.count("o, opt")) {
        cout << _op_desc << endl;
        throw runtime_error("Please specify a valid main program option file.");
    }

    if(_op_vm.count("begin")) {
        has_begin = true;
        begin = _op_vm["begin"].as<double>();
    }

    nframes = _op_vm["nframes"].as<int>();
    first_frame = _op_vm["first-frame"].as<int>();

    CheckInput();

    while (BeginEvaluate(nframes)==true){
        EvaluateFrame();
        nframes--;
    }
    EndEvaluate();
}