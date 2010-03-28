#include "qmapplication.h"

QMApplication::QMApplication()
{}

QMApplication::~QMApplication()
{}

void QMApplication::ParseCommandLine(int argc, char **argv)
{
    namespace po = boost::program_options;

    /// define standard program options
    _op_desc.add_options()
    ("help", "  produce this help message")
    ("crg", boost::program_options::value<string>()->default_value("list_charges.xml"), "  charge unit definitions")
    ("opt", boost::program_options::value<string>()->default_value("main.xml"), "  main program options")
    ("out", boost::program_options::value<string>()->default_value("stateOut.dat"), "  write new state file with this name")
    ("in", boost::program_options::value<string>()->default_value("stateIn.dat"), "  read state file with this name")
    ("nnnames", boost::program_options::value<string>()->default_value("*"), "  List of strings that the concatenation of the two molnames must match to be analyzed")
    ("first-frame", boost::program_options::value<int>()->default_value(0), "  start with this frame")
    ("nframes", boost::program_options::value<int>()->default_value(1), "  process so many frames")
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
        int first_frame = _op_vm["first-frame"].as<int>(); /// starting frame
        int nframes = _op_vm["nframes"].as<int>(); /// number of frames to be processed

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
        
        cout << "Loading qmtopology via state saver." << endl;
        string statefile = _op_vm["in"].as<string>();
        StateSaver loader(_qmtop, statefile,'r');
        string stateout=_op_vm["out"].as<string>();
        StateSaver saver(_qmtop, stateout, 'w');

        loader.Seek(first_frame);
        for (int i=0;i<nframes;i++){
            loader.Load();
            _qmtop.setStep(i);// this is not read by the loader - nor written by the saver (should change)
            EvaluateFrame();
            if (_op_vm.count("out")){
                saver.Save();
            }
        }
        loader.Close();
        saver.Close();
        EndEvaluate();

    }
    catch (std::exception &error) {
        cerr << "an error occured:\n" << error.what() << endl;
    }
}

void QMApplication::HelpText()
{
    //votca::md2qm::HelpTextHeader("unknown program name");
    cout << "no help text available\n\n";
    cout << _op_desc << endl;
    cout << _op_desc_specific << endl;
}

void QMApplication::PrintNbs(string filename){
    ofstream out_nbl;
    out_nbl.open(filename.c_str());
    QMNBList &nblist = _qmtop.nblist();
    if(out_nbl!=0){
        out_nbl << "Neighbours, J(0), J_eff, rate, r_ij, abs(r_ij) [nm]" << endl;
        QMNBList::iterator iter;
        for ( iter  = nblist.begin(); iter != nblist.end() ; ++iter){
            ///Hack!
            if ((*iter)->Js().size() > 0 ){
                out_nbl << "(" << (*iter)->first->getId() << "," << (*iter)->second->getId() << "): ";
                out_nbl << (*iter)->Js()[0] << " " << sqrt((*iter)->calcJeff2()) << " " << (*iter)->rate12() << " ";
                out_nbl << (*iter)->r().getX() << " " << (*iter)->r().getY() << " " << (*iter)->r().getZ() << " ";
                out_nbl << " " << (*iter)->dist() << endl;
            }
        }
    }
    out_nbl.close();
}