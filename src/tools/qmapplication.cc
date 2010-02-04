#include "qmapplication.h"

QMApplication::QMApplication()
{}

QMApplication::~QMApplication()
{}

void QMApplication::setNNnames(string  nnnames){
    Tokenizer tok(nnnames, " ;");
    Tokenizer::iterator itok = tok.begin();
    for (; itok!= tok.end(); ++itok){
        _nnnames.push_back(*itok);
    }
}

bool QMApplication::MatchNNnames(CrgUnit *crg1, CrgUnit* crg2){
    vector <string>::iterator its = _nnnames.begin();
    string namecrg = crg1->getType()->GetName()+string(":")+crg2->getType()->GetName();
    for ( ; its!= _nnnames.end(); ++its){
        if(wildcmp(its->c_str(), namecrg.c_str()) ){
            return true;
        }
    }
    return false;
}

void QMApplication::ParseCommandLine(int argc, char **argv)
{
    namespace po = boost::program_options;

    /// define standard program options
    _op_desc.add_options()
    ("help", "  produce this help message")
    ("crg", boost::program_options::value<string>()->default_value("list_charges.xml"), "  charge unit definitions")
    ("opt", boost::program_options::value<string>()->default_value("main.xml"), "  main program options")
    ("out", boost::program_options::value<string>(), "  write new state file with this name")
    ("nnnames", boost::program_options::value<string>()->default_value("*"), "  List of strings that the concatenation of the two molnames must match to be analyzed")
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

    /// load crg unit definitions from list_charges.xml
    _qmtop.LoadListCharges(_op_vm["crg"].as<string>());

    /// set the nearest neighbor names to be matched
    setNNnames(_op_vm["nnnames"].as<string>());

    /// read in program options from main.xml
    load_property_from_xml(_options, _op_vm["opt"].as<string>());
}

void QMApplication::Run(int argc, char **argv)
{
    try {

        AddSpecificOptions(); /// initialize program-specific parameters
        ParseCommandLine(argc, argv); /// initialize general parameters & read input file
        Initialize(); /// initialize program-specific parameters

        bool has_begin = false; /// was a starting time specified?
        double begin; /// starting time
        int first_frame; /// starting frame
        int nframes = 1; /// number of frames to be processed

        if (_op_vm.count("help")) {
            HelpText();
            return;
        }

        if (!_op_vm.count("opt")) {
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
        StateSaver saver(_qmtop);
        cout << "Loading qmtopology via state saver." << endl;
        string statefile = "falks.dat";
        saver.Load(statefile);

        EvaluateFrame();

        if (_op_vm.count("out")){
           saver.Save(_op_vm["out"].as<string > ());
        }
        else{
            saver.Save(statefile, false);
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
