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
#include "imcio.h"

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
    string name_in, name_out, name_unpack;
    // program options
    po::options_description desc("Allowed options");

    desc.add_options()
      ("in", po::value<string>(&name_in), "files to read")
      ("out", po::value<string>(&name_out), "files to write")
      ("unpack", po::value<string>(&name_unpack), "extract all tables from this file")
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

    ub::vector<double> r;
    ub::vector<double> dS;
    ub::symmetric_matrix<double> gmc;
    vector<string> names;
    vector<RangeParser> ranges;

    imcio_read_index(name_in + ".idx", names, ranges);

    if(vm.count("unpack")) {
        RangeParser *cur_rp;
        Table tbl_in;
        tbl_in.Load(name_unpack);
        vector<string>::iterator iter_name = names.begin();
        vector<RangeParser>::iterator iter_range = ranges.begin();

        while(iter_name != names.end()) {
            cur_rp = &(*iter_range);
            Table tbl;
            for(RangeParser::iterator ir=cur_rp->begin(); ir!=cur_rp->end(); ++ir) {
                tbl.push_back(tbl_in.x(*ir-1), tbl_in.y(*ir-1), 'i');
            }
            tbl.Save(*iter_name + ".dpot.imc");
            ++iter_name;
            ++iter_range;
        }
    } else {
        check_option(desc, vm, "out");
        RangeParser *cur_rp;

        vector<string>::iterator iter_name = names.begin();
        vector<RangeParser>::iterator iter_range = ranges.begin();
        int beg=1;
        int end=0;
        list<int> list;
        
        imcio_read_dS(name_in + ".imc", r, dS);
        imcio_read_matrix(name_in + ".gmc", gmc);

        while(iter_name != names.end()) {
            cur_rp = &(*iter_range);
            for(RangeParser::iterator ir=cur_rp->begin(); ir!=cur_rp->end(); ++ir) {
                for(int i=0; i<gmc.size1(); ++i)
                    if(fabs(gmc(i,*ir-1)) > 1e-8) {
                        list.push_back(*ir-1);
                        end++;
                        break;
                    }                
            }
            RangeParser new_rp;
            new_rp.Add(beg, end);
            beg=end+1;
            *iter_range = new_rp;
            ++iter_name;
            ++iter_range;
        }
        imcio_write_dS(name_out + ".imc", r, dS, &list);
        imcio_write_matrix(name_out + ".gmc", gmc, &list);
        imcio_write_index(name_out + ".idx", names, ranges);
    }

}

