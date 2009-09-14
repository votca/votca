/* 
 * File:   csg_imcrepack.cpp
 * Author: ruehle
 *
 * Created on September 10, 2009, 5:14 PM
 *
 *  csg_imcrepack repacks the matri for the imc update. It either removes zero line entries (pack)
 *  or extracts the dU from the big file (unpack)
 *
 */

#include <tools/table.h>
#include <tools/tokenizer.h>
#include <boost/program_options.hpp>
#include <iostream>

using namespace std;
namespace po = boost::program_options;

void check_option(po::options_description &desc, po::variables_map &vm, const string &option)
{
    if(!vm.count(option)) {
        cout << "csg_imcrepack \n\n";
        cout << desc << endl << "parameter " << option << " is not specified\n";
        exit(1);
    }
}

int main(int argc, char** argv)
{
    string name_in, name_out;
    // program options
    po::options_description desc("Allowed options");

    desc.add_options()
      ("in", po::value<string>(&name_in), "files to read")
      ("out", po::value<string>(&name_out), "files to write")
      ("unpack", "extract all tables")
      ("help", "display help message");

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    }
    catch(po::error err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }

    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_imcrepack \n\n";
        cout << desc << endl;
        return 0;
    }

    check_option(desc, vm, "in");

}

