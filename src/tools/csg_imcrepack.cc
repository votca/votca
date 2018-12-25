/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

/* 
 *  csg_imcrepack repacks the matri for the imc update. It either removes zero line entries (pack)
 *  or extracts the dU from the big file (unpack)
 *
 */

#include <votca/tools/table.h>
#include <votca/tools/tokenizer.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <votca/csg/imcio.h>
#include <votca/csg/version.h>

using namespace votca::csg;
using namespace votca::tools;

void help_text(void)
{
    votca::csg::HelpTextHeader("csg_imcrepack");
    cout << "This program is internally called by inversion scripts to kick out\n"
            "zero entries in matrix for inverse Monte Carlo. It also extracts the\n"
            "single potential updates out of the full solution.\n\n";
}

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
    catch(po::error& err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }

    // does the user want help?
    if (vm.count("help")) {
        help_text();
        cout << desc << endl;
        return 0;
    }

    check_option(desc, vm, "in");

    Eigen::VectorXd r;
    Eigen::VectorXd dS;
    Eigen::MatrixXd gmc;
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
                for(long int i=0; i<gmc.rows(); ++i)
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

