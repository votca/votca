/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <votca/csg/version.h>
#include <iostream>
#include <fstream>


using namespace std;
namespace po = boost::program_options;
using namespace votca::csg;
using namespace votca::tools;

/**
 *  *** Analyze particle distribution as a function of a coordinate ***
 *
 *  This program reads a topology and (set of) trajectory(ies). For every
 *  binned value of a chose coordinate, it outputs the time-averaged number of
 *  particles, listed by particle types.
 */
void help_text(void)
{
    votca::csg::HelpTextHeader("csg_part_dist");
    cout << "This program reads a topology and (set of) trajectory(ies). For every\n"
	"binned value of a chosen coordinate, it outputs the time-averaged number of\n"
	"particles, listed by particle types.\n\n";
}

void check_option(po::options_description &desc, po::variables_map &vm, const string &option)
{
    if(!vm.count(option)) {
        cout << "csg_part_dist \n\n";
        cout << desc << endl << "parameter " << option << " is not specified\n";
        exit(1);
    }
}

int main(int argc, char** argv)
{
    string top_file, trj_file, ptypes_file, out_file, grid, 
	comment, string_tmp, coordinate="z";
    double min, max, step, coord, com(0.);
    vector<int> ptypes;
    ifstream fl_ptypes;
    ofstream fl_out;
    Topology top;
    TopologyReader *reader;
    TrajectoryReader *trajreader;
    int part_type, n_bins, first_frame(0), last_frame(-1), 
	flag_found(0), **p_occ = NULL, n_part(0), frame_id(0),
	analyzed_frames(0);
    bool moreframes(1), not_the_last(1);

    MoleculeContainer::iterator mol;    
    // Load topology+trajectory formats
    TopologyReader::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();

    // read program options
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    
    // let cg_engine add some program options
    desc.add_options()
	("top", po::value<string>(&top_file), "topology file")
	("trj", po::value<string>(&trj_file), "trajectory file")
	("grid", po::value<string>(&grid), "output grid spacing (min:step:max)")
	("out", po::value<string>(&out_file), 
	 "output particle distribution table")
	("ptypes", po::value<string>(&ptypes_file),
	 "particle types to include in the analysis\n"
	 " arg: file - particle types separated by space"
	 "\n default: all particle types")
	("first_frame", po::value<int>(&first_frame), "first frame considered for analysis")
	("last_frame", po::value<int>(&last_frame), "last frame considered for analysis")
	("coord", po::value<string>(&coordinate), "coordinate analyzed ('x', 'y', or 'z' (default))")
	("shift_com", "shift center of mass to zero")
	("comment", po::value<string>(&comment), "store a comment in the output table")
	("help", "produce this help message");
    
    // now read in the command line
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

    check_option(desc, vm, "top");
    check_option(desc, vm, "trj");
    check_option(desc, vm, "out");
    check_option(desc, vm, "grid");

    if (coordinate.compare("x") != 0 && coordinate.compare("y") != 0 && coordinate.compare("z") != 0) {
	cout << "Bad format for coordinate: " << coordinate << endl;
	return 1;
    }

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

    try {
	// Load topology
        reader = TopReaderFactory().Create(vm["top"].as<string>());
        if(reader == NULL) 
            throw std::runtime_error("input format not supported: " + vm["top"].as<string>());
        
        reader->ReadTopology(vm["top"].as<string>(), top);


	// Read the particle types file and save to variable ptypes
	if(vm.count("ptypes")) {
	    fl_ptypes.open(vm["ptypes"].as<string>().c_str());
	    if(!fl_ptypes.is_open())
		throw std::runtime_error("can't open " + vm["ptypes"].as<string>());
	    // build the list of particle types
	    if(fl_ptypes.eof())
		throw std::runtime_error("file " + vm["ptypes"].as<string>() + " is empty");
	    while (!fl_ptypes.eof()) {
		// Not very elegant, but makes sure we don't count the same element twice
		string_tmp = "__";
		fl_ptypes >> string_tmp;
		if (string_tmp != "" && string_tmp != "__") {
		    ptypes.push_back(atoi(string_tmp.c_str()));
		}
	    }
	    fl_ptypes.close();
	    
	}
	else {
	    // Include all particle types
    for(mol=top.Molecules().begin(); mol!=top.Molecules().end();++mol) {
      for(int i=0; i<(*mol)->BeadCount(); ++i) {
        flag_found = 0;
        part_type = atoi((*mol)->getBead(i)->getBeadTypeName().c_str());
        for (size_t j=0; j < ptypes.size(); ++j) {
          if (part_type == ptypes[j]) 
            flag_found = 1;
        }
        if (!flag_found)
          ptypes.push_back(part_type);
      }
    }
	}

	// Allocate array used to store particle occupancy p_occ
	p_occ = (int**) calloc(ptypes.size(),sizeof(int*));
	for (size_t i=0; i < ptypes.size(); ++i)
	    p_occ[i] = (int*) calloc(n_bins,sizeof(int));

	// If we need to shift the center of mass, calculate the number of
	// particles (only the ones that belong to the particle type index
	// ptypes)
  if (vm.count("shift_com")) {
    for(mol=top.Molecules().begin(); mol!=top.Molecules().end();++mol) {
      for(int i=0; i<(*mol)->BeadCount(); ++i) {
        part_type = atoi((*mol)->getBead(i)->getBeadTypeName().c_str());
        for (size_t j=0; j < ptypes.size(); ++j) 
          if (part_type == ptypes[j])
            ++n_part;
      }
    }
  }


	// Now load trajectory
	trajreader = TrjReaderFactory().Create(vm["trj"].as<string>());
	if(trajreader == NULL) 
            throw std::runtime_error("input format not supported: " + vm["trj"].as<string>());
	trajreader->Open(vm["trj"].as<string>());

	// Read the trajectory. Analyze each frame to obtain
	// particle occupancy as a function of coordinate z.	
	while (moreframes) {
	    // Read frame
	    if (frame_id == 0) 
		moreframes = trajreader->FirstFrame(top);
	    else
		moreframes = trajreader->NextFrame(top);

	    // Was this the last frame we read?
	    if (last_frame == -1)
		not_the_last = 1;
	    else {
		if (frame_id <= last_frame)
		    not_the_last = 1;
		else
		    not_the_last = 0;
	    }


	    // Calculate new center of mass position in the direction of 'coordinate'
	    com = 0.;
	    if (vm.count("shift_com")) {
        for(mol=top.Molecules().begin(); mol!=top.Molecules().end();++mol) {
          for(int i=0; i<(*mol)->BeadCount(); ++i) {
            part_type = atoi((*mol)->getBead(i)->getBeadTypeName().c_str());
            for (size_t j=0; j < ptypes.size(); ++j) {
              if (part_type == ptypes[j]) {
                if (coordinate.compare("x") == 0) {
                  com += (*mol)->getBead(i)->getPos().getX();
                }
                else if (coordinate.compare("y") == 0) {
                  com += (*mol)->getBead(i)->getPos().getY();
                }
                else {
                  com += (*mol)->getBead(i)->getPos().getZ();
                }
              }
            }
          }
        }
        com /= n_part;
	    }

	    // Analyze frame
	    if (moreframes && frame_id >= first_frame && not_the_last) {
		++analyzed_frames;
		// Loop over each atom property
    for(mol=top.Molecules().begin(); mol!=top.Molecules().end();++mol) {
      for(int i=0; i<(*mol)->BeadCount(); ++i) {
        part_type = atoi((*mol)->getBead(i)->getBeadTypeName().c_str());
        for (size_t j=0; j < ptypes.size(); ++j) {	
          if (part_type == ptypes[j]) {
            if (coordinate.compare("x") == 0)
              coord = (*mol)->getBead(i)->getPos().getX();
            else if (coordinate.compare("y") == 0)
              coord = (*mol)->getBead(i)->getPos().getY();
            else 
              coord = (*mol)->getBead(i)->getPos().getZ();

            if (coord-com > min && coord-com < max)
              ++p_occ[j][(int)floor((coord-com-min)/step)];
          }
        }
      }
    }
      }
      ++frame_id;
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
	for (size_t i=0; i < ptypes.size(); ++i)
	    fl_out << "type " << ptypes[i] << "\t" << flush;
	fl_out << endl;
	for (int k=0; k < n_bins; ++k) {
	    fl_out << min+k*step << "\t" << flush;
	    for (size_t j=0; j < ptypes.size(); ++j) {
		if (p_occ[j][k] == 0)
		    fl_out << 0 << "\t" << flush;
		else
		    fl_out << p_occ[j][k]/(1.*analyzed_frames) << "\t" << flush;
	    }
	    fl_out << endl;
	}
	fl_out.close();
    }
    catch(std::exception &error) {
	cerr << "An error occured!" << endl << error.what() << endl;
    }

    for (size_t i=0; i < ptypes.size(); ++i)
	free(p_occ[i]);
    free(p_occ);

    cout << "The table was written to " << vm["out"].as<string>() << endl;
    
    return 0;
}

