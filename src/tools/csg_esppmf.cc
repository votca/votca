/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#include <topologyreader.h>
#include <trajectoryreader.h>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include "version.h"
#include <iostream>
#include <fstream>


using namespace std;
namespace po = boost::program_options;

/**
 *  *** Analyze PMF from ESPResSo simulation ***
 *
 *  This program reads an Espresso blockfile, as well as the particle types to
 *  be analyzed. It outputs the binned average particle occupancy (as a function
 *  of the coordinate z), sorted by particle types.
 */
void help_text(void)
{
    votca::csg::HelpTextHeader("csg_esppmf");
    cout << "This program reads an Espresso blockfile, as well as the particle types to\n"
	"be analyzed. It outputs the binned average particle occupancy (as a function\n"
	"of the coordinate z), sorted by particle types.\n\n";
}

void check_option(po::options_description &desc, po::variables_map &vm, const string &option)
{
    if(!vm.count(option)) {
        cout << "csg_esppmf \n\n";
        cout << desc << endl << "parameter " << option << " is not specified\n";
        exit(1);
    }
}

int main(int argc, char** argv)
{
    string in_file, ptypes_file, out_file, grid, comment, string_tmp;
    double min, max, step, **p_occ = NULL, z_coord, sum_weights;
    vector<int> ptypes;
    ifstream fl_ptypes;
    ofstream fl_out;
    Topology top;
    TopologyReader *reader;
    TrajectoryReader *trajreader;
    int part_type, n_bins;
    

    // Load topology+trajectory formats
    TopologyReader::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();

    // read program options
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    
    // let cg_engine add some program options
    desc.add_options()
	("in", po::value<string>(&in_file), 
	 "input ESPResSo topology+trajectory blockfile (*.esp)")
	("ptypes", po::value<string>(&ptypes_file),
	 "particle types to include in the analysis")
	("out", po::value<string>(&out_file), 
	 "output particle type occupancy table")
	("grid", po::value<string>(&grid), "grid spacing (min:step:max)")
	("comment", po::value<string>(&comment), "store a comment in the output table")
	("help", "produce this help message");
    
    // now read in the command line
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
        help_text();
        cout << desc << endl;
        return 0;
    }

    check_option(desc, vm, "in");
    check_option(desc, vm, "ptypes");
    check_option(desc, vm, "out");
    check_option(desc, vm, "grid");

    // Parse grid
    Tokenizer tok(grid, ":");
    vector<string> toks;
    tok.ToVector(toks);
    if(toks.size()!=3) {
	cout << "Wrong range format, use min:step:max\n";
	return 1;
    }
    min = boost::lexical_cast<double>(toks[0]);
    step = boost::lexical_cast<double>(toks[1]);
    max = boost::lexical_cast<double>(toks[2]);
    // Calculate number of bins
    n_bins = (int)((max-min)/(1.*step)+1);

    // Read the particle types file and save to variable ptypes
    try {
	fl_ptypes.open(vm["ptypes"].as<string>().c_str());
	if(!fl_ptypes.is_open())
	    throw std::runtime_error("can't open " + vm["ptypes"].as<string>());
	// build the list of particle types
	if(fl_ptypes.eof())
	    throw std::runtime_error("file " + vm["ptypes"].as<string>() + " is empty");
	while (!fl_ptypes.eof()) {
	    string_tmp = "__";
	    fl_ptypes >> string_tmp;
	    if (string_tmp != "" && string_tmp != "__") {
		ptypes.push_back(atoi(string_tmp.c_str()));
	    }
	}
	fl_ptypes.close();

	// Allocate array used to store particle occupancy p_occ
	p_occ = (double**) calloc(ptypes.size(),sizeof(double));
	for (int i=0; i < ptypes.size(); ++i)
	    p_occ[i] = (double*) calloc(n_bins,sizeof(double));
    }
    catch (std::exception &error) {
	cerr << "An error occured!" << endl << error.what() << endl;
    }



    // try to read the trajectory+topology. Analyze each frame to obtain
    // particle occupancy as a function of coordinate z.
    try {
	// Load topology
        reader = TopReaderFactory().Create(vm["in"].as<string>());
        if(reader == NULL) 
            throw std::runtime_error("input format not supported: " + vm["in"].as<string>());
        
        reader->ReadTopology(vm["in"].as<string>(), top);

	// Now load trajectory
	trajreader = TrjReaderFactory().Create(vm["in"].as<string>());
	if(trajreader == NULL) 
            throw std::runtime_error("input format not supported: " + vm["in"].as<string>());

	trajreader->Open(vm["in"].as<string>());
	
	bool moreframes=1;
	bool firstframe=1;
	while (moreframes) {
	    if (firstframe) {
		moreframes = trajreader->FirstFrame(top);
		firstframe = 0;
	    }
	    else 
		moreframes = trajreader->NextFrame(top);

	    // Analyze frame
	    if (moreframes) {
		// Loop over each atom property - MOVE TO ANALYZE ALL FRAMES
		MoleculeContainer::iterator mol;
		for(mol=top.Molecules().begin(); mol!=top.Molecules().end();++mol) {
		    for(int i=0; i<(*mol)->BeadCount(); ++i) {
			part_type = atoi((*mol)->getBead(i)->getType()->getName().c_str());
			for (int j=0; j < ptypes.size(); ++j) {
			    if (part_type == ptypes[j]) {
				z_coord = (*mol)->getBead(i)->getPos().getZ();
				if (z_coord > min && z_coord < max)
				    ++p_occ[j][(int)floor((z_coord-min)/step)];
				cout << (*mol)->getBeadName(i) << " " << part_type
				     << " z:" << z_coord << endl;
				
			    }
			}
		    }
		}
	    }
	}
	
	trajreader->Close();				
	
    }
    catch(std::exception &error) {
	cerr << "An error occured!" << endl << error.what() << endl;
    }
    
    
    // Output particle occupancy
    try {
	fl_out.open(vm["out"].as<string>().c_str());
	if(!fl_out.is_open())
	    throw std::runtime_error("can't open " + vm["out"].as<string>());

	fl_out << "#z\t" << flush;
	for (int i=0; i < ptypes.size(); ++i)
	    fl_out << "type " << ptypes[i] << "\t" << flush;
	fl_out << endl;
	for (int k=0; k < n_bins; ++k) {
	    fl_out << (k-min)*step << "\t" << flush;
	    // Measure particle occupancy such that the sum of all fields for
	    // each bin is 1
	    sum_weights=0.;
	    for (int j=0; j < ptypes.size(); ++j) {
		sum_weights += p_occ[j][k];
	    }
	    // In case no particle was found, avoid dividingy by zero
	    if (sum_weights == 0)
		sum_weights = 1;
	    for (int j=0; j < ptypes.size(); ++j) {
		fl_out << p_occ[j][k]/sum_weights << "\t" << flush;
	    }
	    fl_out << endl;
	}
	fl_out.close();
    }
    catch(std::exception &error) {
	cerr << "An error occured!" << endl << error.what() << endl;
    }
    

    
    return 0;
}

