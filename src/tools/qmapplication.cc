#include "qmapplication.h"

void QMApplication::ParseCommandLine(int argc, char **argv)
{
    namespace po = boost::program_options;

    /// define standard program options
    _op_desc.add_options()
    ("help", "  produce this help message")
    ("o, opt", boost::program_options::value<string>()->default_value("main.xml"), "  main program options")
    ("out", boost::program_options::value<string>(), "  write new state file")
    ;

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
    try {

        Initialize(); /// initialize program-specific parameters
        ParseCommandLine(argc, argv); /// initialize general parameters & read input file

        bool has_begin = false; /// was a starting time specified?
        double begin; /// starting time
        int first_frame; /// starting frame
        int nframes = 1; /// number of frames to be processed

        if (!_op_vm.count("help")) {
            HelpText();
            return;
        }

        if (!_op_vm.count("o, opt")) {
            cout << _op_desc << endl;
            throw runtime_error("Please specify a valid main program option file.");
        }

        if (_op_vm.count("begin")) {
            has_begin = true;
            begin = _op_vm["begin"].as<double>();
        }

        CheckInput();

        if (!BeginEvaluate()) return;

        /// load qmtop from state saver
        StateSaver saver(*_qmtop);
        saver.Load("state.dat");

        EvaluateFrame();

        if (_op_vm.count("out")){
           saver.Save(_op_vm["out"].as<string > ());
        }
        else{
            saver.Save("state.dat", false);
        }

        EndEvaluate();

    } catch (std::exception &error) {
        cerr << "an error occured:\n" << error.what() << endl;
    }
}

void QMApplication::HelpText()
{
    //votca::md2qm::HelpTextHeader("unknown program name");
    cout << "no help text available\n\n";
    cout << _op_desc << endl;
}
