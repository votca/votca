/* 
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <chrono>
#include <random>
#include <string>

// Boost includes
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

// VOTCA includes
#include <votca/tools/cubicspline.h>
#include <votca/tools/table.h>
#include <votca/tools/linalg.h>

// Local VOTCA includes
#include "votca/csg/nblistgrid.h"
#include "votca/csg/nblistgrid_3body.h"
#include "votca/csg/beadlist.h"

// Local private includes
#include "csg_ml.h"

int main(int argc, char** argv)
{
    CGMachineLearning app;
    return app.Exec(argc, argv);
}


void CGMachineLearning::Initialize(void)
{
    CsgApplication::Initialize();
    AddProgramOptions()
        ("options", boost::program_options::value<string>(), "  options file for machine learning")
        ("trj-force", boost::program_options::value<string>(), "  coarse-grained trajectory containing forces of already known interactions");
}

bool CGMachineLearning::EvaluateOptions()
{
    CsgApplication::EvaluateOptions();
    CheckRequired("trj", "no trajectory file specified");
    CheckRequired("options", "need to specify options file");
    LoadOptions(OptionsMap()["options"].as<string>());

    has_existing_forces_ = false;
    if(OptionsMap().count("trj-force"))
      has_existing_forces_ = true;
    return true;
}

void CGMachineLearning::BeginEvaluate(Topology *top, Topology *)
{
    // Number of CG beads in topology
    _nbeads = top->BeadCount();
    // Set frame counter to zero
    _frame_counter = 0;

    // accuracy for evaluating the difference in bead positions (default 1e-5)
    _dist = 1e-5;
    if (_options.exists("cg.ml.dist")) {
        _dist = _options.get("cg.ml.dist").as<double>();
    }

    // read _nframes from input file
    _nframes = _options.get("cg.ml.frames").as<int>();

    // default setting is _nbeads_per_frame is equal to _nbeads
    _nbeads_per_frame = _nbeads;
    // read _beads_per_nbeads_per_frame from input file
    if (_options.exists("cg.ml.nbeads_per_frame")) {
        _nbeads_per_frame = _options.get("cg.ml.nbeads_per_frame").as<int>();
    }
    // check if _nbeads_per_frame is les or equal read in bead number
    // sanity check
    if ( (_nbeads_per_frame < 0) ||  (_nbeads_per_frame > _nbeads) ) {
        cerr << "\nERROR in CGMachineLearning::BeginEvaluate - _nbeads_per_frame has incorrect value" << endl;
        cerr << "_nbeads_per_frame has to be greater than zero and equal or less the number of CG beads" << endl;
        cerr << "Check your input!" << endl;
        exit(-1);
    }

    // csg_ml is used in training mode, if not specifically given otherwise
    _train = true;
    if (_options.exists("cg.ml.train")) {
        _train = _options.get("cg.ml.train").as<bool>();
    }
    //beads are selected randomly, if not specified otherwise
    _random_selection = true;
    if (_options.exists("cg.ml.random_selection")) {
        _random_selection  = _options.get("cg.ml.random_selection").as<bool>();
    }

    //default regularization parameter
    _lambda = 0.0000000001;
    //read in _lambda
    if (_options.exists("cg.ml.lambda")) {
        _lambda = _options.get("cg.ml.lambda").as<double>();
    }
    std::cout << "\nYou are using VOTCA!\n";
    std::cout << "\nHey, somebody wants to machine learn!\n";
    std::cout << std::endl;

    // initializing non-bonded interactions
    for (votca::tools::Property *prop : _nonbonded) {
        // add mlinfo object to container
        _mls.emplace_back(_mls.size(), prop);    
    }

    //initialize _bead_to_row
    _bead_to_row=Eigen::VectorXi::Zero(_nbeads_per_frame);

    if (_train){
        //initialize _b (_x, _L_K, and _L_T_L_K have to be dynamically adjusted)
        _b=Eigen::VectorXd::Zero(_nbeads_per_frame*_nframes*3);
    }
    else{
        _b=Eigen::VectorXd::Zero(_nbeads_per_frame*3);
        deSerializeMLObjects();
    }

  if (has_existing_forces_) {
    top_force_.CopyTopologyData(top);
    trjreader_force_ =
        TrjReaderFactory().Create(OptionsMap()["trj-force"].as<string>());
    if (trjreader_force_ == nullptr) {
      throw runtime_error(string("input format not supported: ") +
                          OptionsMap()["trj-force"].as<string>());
    }
    // open the trajectory
    trjreader_force_->Open(OptionsMap()["trj-force"].as<string>());
    // read in first frame
    trjreader_force_->FirstFrame(top_force_);
  }
}

CGMachineLearning::MLInfo::MLInfo(int index, votca::tools::Property *options)
{
    // initialize standard data
    MLIndex = index;
    _options = options;
    MLName = options->get("name").value();

    // check if non-bonded 3-body interaction or not (default is no)
    threebody = false;
    if (options->exists("threebody")) {
        threebody = options->get("threebody").as<bool>();
    }

    MLObject = new ML();
    MLObject->setThreebody(threebody);
    
    // initialize threebody_symmetric with true. It is only consider when it is a threebody interaction
    threebody_symmetric = true;

    // check if threebody interaction or not
    if (threebody) {
        std::cout << "Settings for interaction " << MLName << ": " << std::endl;
        std::cout << "You are parametrizing a threebody interaction!" << std::endl << std::endl;
        type1 = options->get("type1").value();
        type2 = options->get("type2").value();
        type3 = options->get("type3").value();
        Eigen::Vector3d sigma;
        sigma(0)=options->get("ml.sigma1").as<double>();
        sigma(1)=options->get("ml.sigma2").as<double>();
        sigma(2)=options->get("ml.sigma3").as<double>();
        MLObject->setSigma(sigma);
        // if type2 is equal to type3, the binning and output tables can be written in a symmetric form to save memory
        if(type2 == type3){
            threebody_symmetric = true;
        }
        // if type2 is not equal to type3, the binning and output tables can not be written in a symmetric form
        if(type2 != type3){
            threebody_symmetric = false;
        }
    } else {
        std::cout << "Settings for interaction " << MLName << ": " << std::endl;
        std::cout << std::endl;
        type1 = options->get("type1").value();
        type2 = options->get("type2").value();
        //initialize with dummy value
        type3 = "dummytype";
        Eigen::Vector3d sigma;
        sigma(0)=options->get("ml.sigma").as<double>();
        //initialize with dummy values
        sigma(1)=1.0;
        sigma(2)=1.0;
        MLObject->setSigma(sigma);
    }

    // read in the values for output
    min_out = options->get("ml.min").as<double>();
    max_out = options->get("ml.max").as<double>();
    std::cout << "Cutoff distance for the ML parametrization: " << max_out << std::endl << std::endl;
    // initialize grid size for sample predictions
    //default value for output grid size
    dx_out = 0.05;
    if (options->exists("ml.out_step")) {
        dx_out = options->get("ml.out_step").as<double>();
    }
    // number of output grid points
    num_out = 1 + (int)((max_out - min_out) / dx_out);
    // check if threebody interaction or not
    if (threebody) {
        std::cout << "Grid parameters for force predictions on sample triplets: " << std::endl << std::endl;
        std::cout << "Distances r_12 and r_13: min_out " << min_out << ", max_out: " << max_out << ", dx_out: " << dx_out << ", num_out: " << num_out << std::endl;
        std::cout << "Angles: min_out 0.0, max_out: 180.0, dtheta_out: " << 180/(num_out-1) << ", num_out: " << num_out << std::endl << std::endl;
    } else {
        std::cout << "Grid parameters for output table of pair force: " << std::endl << std::endl;
        std::cout << "Distances r: min_out " << min_out << ", max_out: " << max_out << ", dx_out: " << dx_out << ", num_out: " << num_out << std::endl << std::endl;
    }
    //default value of d is 0.0 which means no switching function is applied
    d = 0.0;
    if (options->exists("ml.d")) {
        d = options->get("ml.d").as<double>();
        std::cout << "Switching function applied with width of transition region d: " << d << std::endl;
    }
    else{
        std::cout << "No switching function applied" << std::endl;
    }
    std::cout << std::endl;

    // initialize sum of test errors
    test_error_sum = 0.0;

    //to be changed
    //check if threebody interaction or not
    if (threebody) {
        result=Eigen::MatrixXd::Zero(5*3,num_out*3);
        resulttheta=Eigen::MatrixXd::Zero(4*3,num_out*3);
        error=Eigen::MatrixXd::Zero(5*3,num_out*3);
        errortheta=Eigen::MatrixXd::Zero(4*3,num_out*3);
    } else {
        result=Eigen::MatrixXd::Zero(1,num_out*3);
        error=Eigen::MatrixXd::Zero(1,num_out*3);
        resulttheta=Eigen::MatrixXd::Zero(0,0);
        errortheta=Eigen::MatrixXd::Zero(0,0);
    }

    //true if ml output table should be generated
    output_table = false;
    if (options->exists("ml.output_table")) {
        output_table = options->get("ml.output_table").as<bool>();
    }
    //the default number of grid points for the output table
    num_table = 20;
    if (options->exists("ml.N_table")) {
        num_table = options->get("ml.N_table").as<int>();
    }
    // initialize grid size for force output table
    dx_table = ((max_out - min_out) / (num_table));
    // option only relevant if threebody interaction
    if ( (output_table == true) && (threebody) ){
        //if it is a symmetric threebody interaction, less table entries are required
        if (threebody_symmetric == true){
            // only then resize Eigen vectors to store output table entries
            // number of lines in output table is num_table*num_table*(num_table+1) as number of sampled theta values is 2*num_table
            output1_1=Eigen::VectorXd::Zero(num_table*num_table*(num_table+1));
            output1_2=Eigen::VectorXd::Zero(num_table*num_table*(num_table+1));
            output2_1=Eigen::VectorXd::Zero(num_table*num_table*(num_table+1));
            output2_2=Eigen::VectorXd::Zero(num_table*num_table*(num_table+1));
            output3_1=Eigen::VectorXd::Zero(num_table*num_table*(num_table+1));
            output3_2=Eigen::VectorXd::Zero(num_table*num_table*(num_table+1));
        }
        //if it is not a symmetric threebody interaction, more table entries are required
        if (threebody_symmetric == false){
            // only then resize Eigen vectors to store output table entries
            // number of lines in output table is 2*num_table*num_table*num_table as number of sampled theta values is 2*num_table
            output1_1=Eigen::VectorXd::Zero(2*num_table*num_table*num_table);
            output1_2=Eigen::VectorXd::Zero(2*num_table*num_table*num_table);
            output2_1=Eigen::VectorXd::Zero(2*num_table*num_table*num_table);
            output2_2=Eigen::VectorXd::Zero(2*num_table*num_table*num_table);
            output3_1=Eigen::VectorXd::Zero(2*num_table*num_table*num_table);
            output3_2=Eigen::VectorXd::Zero(2*num_table*num_table*num_table);
        }
        std::cout << "You are generating a force output table!" << std::endl << std::endl;
        std::cout << "Grid parameters for the output table: " << std::endl << std::endl;
        std::cout << "Distances r_12 and r_13: min_out " << (min_out+0.5*dx_table) << ", max_out: " << (max_out-0.5*dx_table) << ", dx_table: " << dx_table << ", num_table: " << num_table << std::endl;
        std::cout << "Angles: min_theta: " << (0.0+0.5*(180.0/(num_table*2))) << ", max_theta: " << (180.0-0.5*(180.0/(num_table*2))) << ", dtheta_table: " << (180.0/(num_table*2)) << ", num_theta: " << (num_table*2) << std::endl << std::endl;
    }
    
    //Check, if binning is switched on for this interaction, by default it is switched off
    binning = false;
    if (options->exists("ml.binning")) {
        binning  = options->get("ml.binning").as<bool>();
        std::cout << "You are using the ML code with covariant meshing!" << std::endl << std::endl;
    }
    
    //these variables and vectors are only initialized when binning
    if (binning == true){
        //look if option min_theta exists, otherwise set it to zero
        min_theta = 0.0;
        if (options->exists("ml.min_theta")) {
            min_theta = options->get("ml.min_theta").as<double>();
        }

        //the default number of bins is 20
        num_bins = 20;
        if (options->exists("ml.N_bins")) {
            num_bins = options->get("ml.N_bins").as<int>();
        }

        //number of bins for theta (only relevant if threebody interaction, otherwise ignored)
        num_theta = num_bins*2;
        if (options->exists("ml.N_theta")) {
            num_theta = options->get("ml.N_theta").as<int>();
        }

        //the default number of bins for the second grid is 0 (so far only implemented for threebody interaction)
        num_bins2 = 0;
        if (options->exists("ml.N_bins2")) {
            num_bins2 = options->get("ml.N_bins2").as<int>();
        }

        //look if option min2 exists (minimum pair distance for second grid)
        // (so far only relevant if threebody interaction, otherwise ignored)
        min_out2 = max_out;
        if (options->exists("ml.min2")) {
            min_out2 = options->get("ml.min2").as<double>();
        }

        //if second grid, min_out2 = max_out, if not min_out2 < max_out (so far only implemented for threebody interaction)
        dx_bins = (min_out2-min_out)/num_bins;

        dx_bins2 = 0;
        //if second grid (so far only implemented for threebody interaction)
        if (min_out2 < max_out){
            dx_bins2 = (max_out-min_out2)/num_bins2;
        }

        //default value of smear_scale is 1.8
        smear_scale = 1.8;
        if (options->exists("ml.smear_scale")) {
            smear_scale = options->get("ml.smear_scale").as<double>();
        }

        //The binning is done using Gaussian smearing, if not given otherwise
        gaussian_smearing = true;
        if (options->exists("ml.gaussian_smearing")) {
            gaussian_smearing  = options->get("ml.gaussian_smearing").as<bool>();
        }

        std::cout << "Grid parameters for the covariant meshing: " << std::endl << std::endl;

        //check if threebody interaction or not
        if (threebody) {
            //if only one grid
            if(min_out2 == max_out){
                std::cout << "Distances r_12 and r_13: min_out " << (min_out+0.5*dx_bins) << ", max_out: " << (max_out-0.5*dx_bins) << ", dx_bins: " << dx_bins << ", num_bins: " << num_bins << std::endl;
                std::cout << "Angles: min_theta: " << (min_theta+0.5*(180.0-min_theta)/(num_theta)) << ", max_theta: " << (180.0-0.5*(180.0-min_theta)/(num_theta)) << ", dtheta_out: " << (180.0-min_theta)/(num_theta) << ", num_theta: " << num_theta << std::endl;
                std::cout << "The total number of grid points for the (three dimensional) covariant meshing is: " << num_bins*num_theta*(num_bins+1)/2 << std::endl << std::endl;
            }
            //if second grid
            if(min_out2 < max_out){
                std::cout << "You are using a combination of two grids with two different grid parameters" << std::endl << std::endl;
                std::cout << "Grid parameters for the first grid: " << std::endl << std::endl;
                std::cout << "Distances r_12 and r_13: min_out " << (min_out+0.5*dx_bins) << ", max_out: " << (min_out2-0.5*dx_bins) << ", dx_bins: " << dx_bins << ", num_bins: " << num_bins << std::endl;
                std::cout << "Angles: min_theta: " << (min_theta+0.5*(180.0-min_theta)/(num_theta)) << ", max_theta: " << (180.0-0.5*(180.0-min_theta)/(num_theta)) << ", dtheta_out: " << (180.0-min_theta)/(num_theta) << ", num_theta: " << num_theta << std::endl;
                std::cout << "Grid parameters for the second grid: " << std::endl << std::endl;
                std::cout << "Distances r_12 and r_13: min_out2 " << (min_out2+0.5*dx_bins2) << ", max_out2: " << (max_out-0.5*dx_bins2) << ", dx_bins2: " << dx_bins2 << ", num_bins2: " << num_bins2 << std::endl;
                std::cout << "Angles: min_theta: " << (min_theta+0.5*(180.0-min_theta)/(num_theta)) << ", max_theta: " << (180.0-0.5*(180.0-min_theta)/(num_theta)) << ", dtheta_out: " << (180.0-min_theta)/(num_theta) << ", num_theta: " << num_theta << std::endl;
                std::cout << "The total number of grid points for the (three dimensional) covariant meshing is: " << (num_bins+num_bins2)*num_theta*(num_bins+num_bins2+1)/2 << std::endl << std::endl;
            }
        } else {
                std::cout << "Distances r: min_out " << (min_out+0.5*dx_bins) << ", max_out: " << (max_out-0.5*dx_bins) << ", dx_bins: " << dx_bins << ", num_bins: " << num_bins << std::endl;
                std::cout << "The total number of grid points for the (three dimensional) covariant meshing is: " << num_bins << std::endl << std::endl;
        }

        //if Gaussian smearing is applied
        if (gaussian_smearing == true){
            std::cout << "You are using Gaussian smearing with smear_scale: " << smear_scale << std::endl;
        }

        //check if threebody interaction or not to use correct grid size and initialization
        if (threebody) {
            double out_12,out_13,theta,dtheta;
            int count;
            dtheta = (180.0-min_theta)/(num_theta);

            //if it is a symmetric threebody interaction, less table entries are required for the binning
            if (threebody_symmetric == true){
                binned_structures=Eigen::VectorXd::Zero(num_bins*num_theta*(num_bins+1)/2*9);

                out_12 = min_out+0.5*dx_bins;
                count = 0;

                for (int i=0; i<num_bins; ++i){
                    out_13 = out_12;
                    for (int j=i; j<num_bins; ++j){
                        theta = min_theta+0.5*dtheta;
                        for (int k=0; k<num_theta; ++k){
                            binned_structures[count*9 + 0] = 0.0;
                            binned_structures[count*9 + 1] = 0.0;
                            binned_structures[count*9 + 2] = out_12;
                            binned_structures[count*9 + 3] = 0.0;
                            binned_structures[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                            binned_structures[count*9 + 5] = out_13*cos(theta * M_PI / 180.0);
                            binned_structures[count*9 + 6] = binned_structures[count*9 + 3] - binned_structures[count*9 + 0];
                            binned_structures[count*9 + 7] = binned_structures[count*9 + 4] - binned_structures[count*9 + 1];
                            binned_structures[count*9 + 8] = binned_structures[count*9 + 5] - binned_structures[count*9 + 2];
                            theta += dtheta;
                            ++count;
                        }
                        out_13 += dx_bins;
                    }
                    out_12 += dx_bins;
                }

                //only second grid if min_out2 < max_out
                if (min_out2 < max_out){
                    binned_structures=Eigen::VectorXd::Zero((num_bins+num_bins2)*num_theta*(num_bins+num_bins2+1)/2*9);

                    out_12 = min_out+0.5*dx_bins;
                    count = 0;

                    for (int i=0; i<(num_bins+num_bins2); ++i){
                        out_13 = out_12;
                        for (int j=i; j<(num_bins+num_bins2); ++j){
                            theta = min_theta+0.5*dtheta;
                            for (int k=0; k<num_theta; ++k){
                                binned_structures[count*9 + 0] = 0.0;
                                binned_structures[count*9 + 1] = 0.0;
                                binned_structures[count*9 + 2] = out_12;
                                binned_structures[count*9 + 3] = 0.0;
                                binned_structures[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                                binned_structures[count*9 + 5] = out_13*cos(theta * M_PI / 180.0);
                                binned_structures[count*9 + 6] = binned_structures[count*9 + 3] - binned_structures[count*9 + 0];
                                binned_structures[count*9 + 7] = binned_structures[count*9 + 4] - binned_structures[count*9 + 1];
                                binned_structures[count*9 + 8] = binned_structures[count*9 + 5] - binned_structures[count*9 + 2];
                                theta += dtheta;
                                ++count;
                            }
         	            if ( j<(num_bins-1) ){
                                out_13 += dx_bins;
           	            }
	                    if ( j==(num_bins-1) ){
                                out_13 += 0.5*dx_bins;
                                out_13 += 0.5*dx_bins2;
	                    }
	                    if ( j>(num_bins-1) ){
                                out_13 += dx_bins2;
	                    }
                        }
	                if ( i<(num_bins-1) ){
                            out_12 += dx_bins;
	                }
	                if ( i==(num_bins-1) ){
                            out_12 += 0.5*dx_bins;
                            out_12 += 0.5*dx_bins2;
	                }
	                if ( i>(num_bins-1) ){
                            out_12 += dx_bins2;
	                }
                    }
                }
            }
            //if it is not symmetric threebody interaction, more table entries are required for the binning
            if (threebody_symmetric == false){
                binned_structures=Eigen::VectorXd::Zero(num_bins*num_theta*num_bins*9);

                out_12 = min_out+0.5*dx_bins;
                count = 0;

                for (int i=0; i<num_bins; ++i){
                    out_13 = min_out+0.5*dx_bins;;
                    for (int j=0; j<num_bins; ++j){
                        theta = min_theta+0.5*dtheta;
                        for (int k=0; k<num_theta; ++k){
                            binned_structures[count*9 + 0] = 0.0;
                            binned_structures[count*9 + 1] = 0.0;
                            binned_structures[count*9 + 2] = out_12;
                            binned_structures[count*9 + 3] = 0.0;
                            binned_structures[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                            binned_structures[count*9 + 5] = out_13*cos(theta * M_PI / 180.0);
                            binned_structures[count*9 + 6] = binned_structures[count*9 + 3] - binned_structures[count*9 + 0];
                            binned_structures[count*9 + 7] = binned_structures[count*9 + 4] - binned_structures[count*9 + 1];
                            binned_structures[count*9 + 8] = binned_structures[count*9 + 5] - binned_structures[count*9 + 2];
                            theta += dtheta;
                            ++count;
                        }
                        out_13 += dx_bins;
                    }
                    out_12 += dx_bins;
                }

                //only second grid if min_out2 < max_out
                if (min_out2 < max_out){
                    binned_structures=Eigen::VectorXd::Zero((num_bins+num_bins2)*num_theta*(num_bins+num_bins2)*9);

                    out_12 = min_out+0.5*dx_bins;
                    count = 0;

                    for (int i=0; i<(num_bins+num_bins2); ++i){
                        out_13 = min_out+0.5*dx_bins;
                        for (int j=0; j<(num_bins+num_bins2); ++j){
                            theta = min_theta+0.5*dtheta;
                            for (int k=0; k<num_theta; ++k){
                                binned_structures[count*9 + 0] = 0.0;
                                binned_structures[count*9 + 1] = 0.0;
                                binned_structures[count*9 + 2] = out_12;
                                binned_structures[count*9 + 3] = 0.0;
                                binned_structures[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                                binned_structures[count*9 + 5] = out_13*cos(theta * M_PI / 180.0);
                                binned_structures[count*9 + 6] = binned_structures[count*9 + 3] - binned_structures[count*9 + 0];
                                binned_structures[count*9 + 7] = binned_structures[count*9 + 4] - binned_structures[count*9 + 1];
                                binned_structures[count*9 + 8] = binned_structures[count*9 + 5] - binned_structures[count*9 + 2];
                                theta += dtheta;
                                ++count;
                            }
         	            if ( j<(num_bins-1) ){
                                out_13 += dx_bins;
           	            }
	                    if ( j==(num_bins-1) ){
                                out_13 += 0.5*dx_bins;
                                out_13 += 0.5*dx_bins2;
	                    }
	                    if ( j>(num_bins-1) ){
                                out_13 += dx_bins2;
	                    }
                        }
	                if ( i<(num_bins-1) ){
                            out_12 += dx_bins;
	                }
	                if ( i==(num_bins-1) ){
                            out_12 += 0.5*dx_bins;
                            out_12 += 0.5*dx_bins2;
	                }
	                if ( i>(num_bins-1) ){
                            out_12 += dx_bins2;
	                }
                    }
                }
            }
        } else {
            binned_structures=Eigen::VectorXd::Zero(num_bins*3);
            double out_12;
            out_12=min_out+dx_bins*0.5;
            for (int i=0; i<num_bins; ++i){
                binned_structures[i*3] = 0.0;
                binned_structures[i*3 + 1] = 0.0;
                binned_structures[i*3 + 2] = out_12;
                out_12 += dx_bins;
            }
        }
    }
    //else they are initialized with default values to prevent uninitialized variables
    else{
        min_out2 = max_out;
        min_theta = 0.0;
        num_bins = 20;
        num_theta = num_bins*2;
        num_bins2 = 0;
        dx_bins = (min_out2-min_out)/num_bins;
        dx_bins2 = 0;
        smear_scale = 1.0;
        gaussian_smearing = true;
    }
}

void CGMachineLearning::EndEvaluate()
{
    // sanity check, to be done    
    std::cout << "\nWe are done, thank you very much!" << std::endl;

  if(has_existing_forces_) {
    trjreader_force_->Close();
  }
}

void CGMachineLearning::SerializeMLObjects(){
    string ML_DIR = "MLObjects";
    string file_extension = ".ml";
    string file_name;
    string ml_id;
    //create directory
    boost::filesystem::create_directories(ML_DIR);
    // save ML objects
    for (MLInfo &mlinfo : _mls) {
        ml_id =  boost::lexical_cast<string>(mlinfo.MLIndex);
        file_name = ml_id+file_extension;
        std::ofstream ofs( (ML_DIR + "/" + file_name).c_str() );
        boost::archive::binary_oarchive oa( ofs );
        oa << mlinfo.MLObject;
    }
}

void CGMachineLearning::deSerializeMLObjects(){

    string ML_DIR = "MLObjects";
    string file_extension = ".ml";
    string file_name;
    string ml_id;
    // load ML objects
    for (MLInfo &mlinfo : _mls) {
        ml_id =  boost::lexical_cast<string>(mlinfo.MLIndex);
        file_name = ml_id+file_extension;
        std::ifstream ifs( (ML_DIR + "/" + file_name).c_str() );
        boost::archive::binary_iarchive ia( ifs );
        ia >> mlinfo.MLObject;
        std::cout << "_sigma deserialize: " << mlinfo.MLObject->getSigma() << std::endl;
        std::cout << "_bead_number deserialize: " << mlinfo.MLObject->getBeadNumber() << std::endl;
        std::cout << "_struct_number deserialize: " << mlinfo.MLObject->getStructNumber() << std::endl;
    }
}

void CGMachineLearning::WriteOutFilesTrain()
{
    std::cout << "Write out files" << std::endl;

    string file_name_base;
    string file_extension_1 = ".force";
    string file_extensiontheta_1 = ".forcetheta";
    string file_extension_2;
    string file_extension_3;
    string file_extension_4;
    string file_name;
    votca::tools::Table force_tab;
    double out_x;
    double dtheta;
    double out_theta;

    // table without error column
    force_tab.SetHasYErr(true);

    for (MLInfo &mlinfo : _mls) {
        //check if threebody
        if (mlinfo.threebody) {
            // construct meaningful outfile name
            file_name_base = mlinfo.MLName;
            for (unsigned int i=0; i<5; ++i){
                for (unsigned int j=0; j<3; ++j){
                    for (unsigned int k=0; k<3; ++k){
                        file_extension_2 = std::to_string(i+1);
                        file_extension_3 = std::to_string(j+1);
                        if (k == 0)
                            file_extension_4 = "_x";
                        if (k == 1)
                            file_extension_4 = "_y";
                        if (k == 2)
                            file_extension_4 = "_z";
                        file_name = file_name_base + file_extension_1 + file_extension_2 + "_" + file_extension_3 + file_extension_4;
                        // resize table
                        force_tab.resize(mlinfo.num_out);

                        // print output file names on stdout
                        std::cout << "Updating file: " << file_name << std::endl;

                        // first output point = first grid point
                        out_x = mlinfo.min_out;
                        // loop over output grid
                        for (int l = 0; l < mlinfo.num_out; l++) {
                            // put point, result, flag and accuracy at point out_x into the table
                            force_tab.set(l, out_x, mlinfo.result(3*i+j,l*3+k), 'i', mlinfo.error(3*i+j,l*3+k));
                            // update out_x for the next iteration
                            out_x += mlinfo.dx_out;
                        }
                        // save table in the file
                        force_tab.Save(file_name);
                        // clear the table for the next spline
                        force_tab.clear();
                    }
                }
            }
            dtheta = 180.0/(mlinfo.num_out-1);
            for (unsigned int i=0; i<4; ++i){
                for (unsigned int j=0; j<3; ++j){
                    for (unsigned int k=0; k<3; ++k){
                        file_extension_2 = std::to_string(i+1);
                        file_extension_3 = std::to_string(j+1);
                        if (k == 0)
                            file_extension_4 = "_x";
                        if (k == 1)
                            file_extension_4 = "_y";
                        if (k == 2)
                            file_extension_4 = "_z";
                        file_name = file_name_base + file_extensiontheta_1 + file_extension_2 + "_" + file_extension_3 + file_extension_4;
                        // resize table
                        force_tab.resize(mlinfo.num_out);

                        // print output file names on stdout
                        std::cout << "Updating file: " << file_name << std::endl;

                        // first output point = first grid point
                        out_theta = 0.0;
                        // loop over output grid
                        for (int l = 0; l < mlinfo.num_out; l++) {
                            // put point, result, flag and accuracy at point out_x into the table
                            force_tab.set(l, out_theta, mlinfo.resulttheta(3*i+j,l*3+k), 'i', mlinfo.errortheta(3*i+j,l*3+k));
                            // update out_theta for the next iteration
                            out_theta += dtheta;
                        }
                        // save table in the file
                        force_tab.Save(file_name);
                        // clear the table for the next spline
                        force_tab.clear();
                    }
                }
            }
        } else {
            // construct meaningful outfile name
            file_name_base = mlinfo.MLName;
            file_name = file_name_base + file_extension_1;
            // resize table
            force_tab.resize(mlinfo.num_out);

            // print output file names on stdout
            std::cout << "Updating file: " << file_name << std::endl;

            // first output point = first grid point
            out_x = mlinfo.min_out;
            // loop over output grid
            for (int l = 0; l < mlinfo.num_out; l++) {
                // put point, result, flag and accuracy at point out_x into the table
                force_tab.set(l, out_x, -mlinfo.result(0,l*3), 'i', mlinfo.error(0,l*3));
                // update out_x for the next iteration
                out_x += mlinfo.dx_out;
            }
            // save table in the file
            force_tab.Save(file_name);
            // clear the table for the next spline
            force_tab.clear();
        }
    }
}

void CGMachineLearning::EvalConfiguration(Topology *conf, Topology *)
{
    if(conf->BeadCount() == 0)
        throw std::runtime_error("CG Topology has 0 beads, check your mapping file!");
    if(has_existing_forces_) {
        if(conf->BeadCount() != top_force_.BeadCount())
            throw std::runtime_error("number of beads in topology and force topology does not match");
        for(votca::Index i=0; i<conf->BeadCount(); ++i) {
            conf->getBead(i)->F() -= top_force_.getBead(i)->getF();
            Eigen::Vector3d d = conf->getBead(i)->getPos() - top_force_.getBead(i)->getPos();
            if(d.norm() > _dist){//default is 1e-5, otherwise it can be a too strict criterion
                throw std::runtime_error("One or more bead positions in mapped and reference force trajectory differ by more than 1e-5");
            }
        }
    }

    //if _nbeads_per_frame < _nbeads and _random_selection == false, choose the first _nbeads_per_frame of the correct bead type per interaction
    //at the moment limited to one interaction (first interaction in container _mls)!!!
    if ((_nbeads_per_frame < _nbeads) && (_random_selection == false)){
        //go through the first _nbeads_per_frame atoms
        //counter for beads of that type
        votca::Index count_beads_of_type1 = 0;
        int iatom = 0;
        while ((count_beads_of_type1 < _nbeads_per_frame) && (iatom < _nbeads)){
            if (conf->getBead(iatom)->getType() == _mls.front().type1){
                //if bead is of type1 (first beadtype of first interaction in _mls!!),
                //set _bead_to_row(count_beads_of_type) to  iatom
               _bead_to_row(count_beads_of_type1) = iatom;
               //now increase counter count_beads_of_type1, as now a bead of type1 has been found
               ++count_beads_of_type1;
            }
            //increase iatom counter in each iteration of the loop
            ++iatom;
        }
        //now check, if enough beads of the specific type have been found
        if (count_beads_of_type1 < _nbeads_per_frame){
        throw std::runtime_error("Not enough beads of type1 in the configuration");
        }
    }
    //if _nbeads_per_frame < _nbeads and _random_selection == true, choose randomly _nbeads_per_frame of the correct bead type per interaction
    //at the moment limited to one interaction (first interaction in container _mls)!!!
    if ((_nbeads_per_frame < _nbeads) && (_random_selection == true)){
        // construct a trivial random generator engine from a time-based seed
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator (seed);
        uniform_int_distribution<int> distribution(0,_nbeads-1);
        for (int iatom = 0; iatom < _nbeads_per_frame; ++iatom) {
            //fill _bead_to_row with integer random numbers between 0 and (_nbeads_per_frame-1)
            //check that no bead occurs twice in the list
            //Preliminary: so far this only works correctly for parametrizing one interaction at a time!!!
            //now also check that the bead found by the random search is indeed of type1!
            bool bead_in_list = true;
            bool bead_of_type1 = false;
            int bead_index;
            while ( (bead_in_list == true) || (bead_of_type1 == false) ){
                bead_in_list = false;
                bead_of_type1 = true;
                bead_index = distribution(generator);
                if (conf->getBead(bead_index)->getType() != _mls.front().type1){
                    bead_of_type1 = false;
                }
                if (bead_of_type1 == true){
                    for (int i = 0; i<(iatom+1); ++i){
                        if (_bead_to_row(i) == bead_index){
                            bead_in_list = true;
                        }
                    }
                }
            }
            //now set index of _bead_to_row(iatom)
            _bead_to_row(iatom) = bead_index;
        }
    }
    //otherwise fill _bead_to_row with integer random numbers between 0 and (_nbeads_per_frame-1)
    if(_nbeads_per_frame == _nbeads){
        for (int iatom = 0; iatom < _nbeads_per_frame; ++iatom) {
            _bead_to_row(iatom)=iatom;
        }
    }

    // loop for the force vector:
    // hack, change the Has functions..
    if (conf->getBead(0)->HasF()) {
        Eigen::Vector3d Force;
        //read in _nbeads_per_frame energies
        for (int iatom = 0; iatom < _nbeads_per_frame; ++iatom) {
            //_bead_to_row connects the bead number with the position in the list
            Force = conf->getBead(_bead_to_row(iatom))->getF();
            if (_train){
            	_b((_nbeads_per_frame * _frame_counter + iatom)*3) = Force.x();
            	_b((_nbeads_per_frame * _frame_counter + iatom)*3 + 1) = Force.y();
            	_b((_nbeads_per_frame * _frame_counter + iatom)*3 + 2) = Force.z();
            }
            else{
            	_b(iatom*3) = Force.x();
            	_b(iatom*3 + 1) = Force.y();
            	_b(iatom*3 + 2) = Force.z();
            }
        }
    } else {
        cerr << "\nERROR in csg_ml::EvalConfiguration - No forces in configuration!\n" << endl;
        exit(-1);
    }

    for (MLInfo &mlinfo : _mls) {
        if (_train){
            // non-bonded interaction
            // check if threebody interaction or not
            if (mlinfo.threebody) {
                EvalNonbondedTrain_Threebody(conf, &mlinfo);
            } else {
                EvalNonbondedTrain(conf, &mlinfo);
            }
        } else {
            // non-bonded interaction testing
            // check if threebody interaction or not
            if (mlinfo.threebody) {
                EvalNonbondedTest_Threebody(conf, &mlinfo);
            } else {
                EvalNonbondedTest(conf, &mlinfo);
            }
        }
    }

    // update the frame counter
    _frame_counter += 1;

    if (_frame_counter % _nframes == 0) { // at this point we processed _nframes frames,
        // is csg_ml run in training or test mode?
        if (_train){
            // solve ML equations and accumulate the result
            AccumulateDataTrain();
            SerializeMLObjects();
            // write results to output files
            WriteOutFilesTrain();
        }
        else{
            // combine test results for different frames
            AccumulateDataTest();
        }
        //go through all interactions and write out output tables
        for (MLInfo &mlinfo : _mls) {
            if (mlinfo.output_table){
                //so far only implemented for threebody interactions
                if (mlinfo.threebody) {
                    //calculate table fo that interaction
                    EvaluateTable_Threebody(&mlinfo);
                    //write table for that interaction
                    WriteTable_Threebody(&mlinfo);
                }
            }
        }
    }
    if(has_existing_forces_)
        trjreader_force_->NextFrame(top_force_);
}

void CGMachineLearning::EvalNonbondedTrain_Threebody(Topology *conf, MLInfo *mlinfo)
{
    // generate the neighbour list
    NBList_3Body *nb;

    bool gridsearch=false;

    if(_options.exists("cg.nbsearch")) {
        if(_options.get("cg.nbsearch").as<string>() == "grid")
            gridsearch=true;
        else if(_options.get("cg.nbsearch").as<string>() == "simple")
            gridsearch=false;
        else throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
    if(gridsearch)
        nb = new NBListGrid_3Body();
    else
        nb = new NBList_3Body();

    nb->setCutoff(mlinfo->_options->get("ml.max").as<double>()); // implement different cutoffs for different interactions!

    // generate the bead lists
    BeadList beads1, beads2, beads3;
    beads1.Generate(*conf, mlinfo->type1);
    beads2.Generate(*conf, mlinfo->type2);
    beads3.Generate(*conf, mlinfo->type3);

    // Generate the 3body neighbour lists
    // check if type2 and type3 are the same
    if (mlinfo->type2 == mlinfo->type3) {
        // if then type2 and type1 are the same, all three types are the same
        // use the Generate function for this case
        if (mlinfo->type1 == mlinfo->type2) {
            nb->Generate(beads1, true);
        }
        // else use the Generate function for type2 being equal to type3 (and type1 being different)
        if (mlinfo->type1 != mlinfo->type2) {
            nb->Generate(beads1, beads2, true);
        }
    }
    // If type2 and type3 are not the same, use the Generate function for three different bead types
    // (Even if type1 and type2 or type1 and type3 are the same, the Generate function for two different beadtypes
    // is only applicable for the case that type2 is equal to type3
    if (mlinfo->type2 != mlinfo->type3) {
        nb->Generate(beads1, beads2, beads3, true);
    }

    // iterate over all triples
    for (BeadTriple *triple : *nb){
        int iatom = triple->bead1()->getId();
        int jatom = triple->bead2()->getId();
        int katom = triple->bead3()->getId();
        Eigen::Vector3d var1 = triple->r12();
        Eigen::Vector3d var2 = triple->r13();
        Eigen::Vector3d var3 = triple->r23();
        double d12 = triple->dist12();
        double d13 = triple->dist13();
        //double d23 = triple->dist23();

        //check if iatom of this triple is in list _bead_to_row
        for (int i=0; i<_nbeads_per_frame; ++i){
            if (iatom == _bead_to_row(i)){
                // if type2 is equal to type3 ordering is such that d13 is always greater than d12
                if (mlinfo->threebody_symmetric == true){
                    if (d13 > d12){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,1);
                    }
                    if (d12 > d13){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var2,var1,-var3,1);
                    }
                }
                // if type2 is not equal to type3 d13 can be smaller and greater than d12
                if (mlinfo->threebody_symmetric == false){
                    conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,1);
                }
            }
            if (jatom == _bead_to_row(i)){
                // if type2 is equal to type3 ordering is such that d13 is always greater than d12
                if (mlinfo->threebody_symmetric == true){
                    if (d13 > d12){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,2);
                    }
                    if (d12 > d13){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var2,var1,-var3,3);
                    }
                }
                // if type2 is not equal to type3 d13 can be smaller and greater than d12
                if (mlinfo->threebody_symmetric == false){
                    conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,2);
                }
            }
            if (katom == _bead_to_row(i)){
                // if type2 is equal to type3 ordering is such that d13 is always greater than d12
                if (mlinfo->threebody_symmetric == true){
                    if (d13 > d12){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,3);
                    }
                    if (d12 > d13){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var2,var1,-var3,2);
                    }
                }
                // if type2 is not equal to type3 d13 can be smaller and greater than d12
                if (mlinfo->threebody_symmetric == false){
                    conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,3);
                }
            }
        }
    }

    int nbeads_old,ntriples_old,nbeads_total,ntriples_total;

    //check, if binning is applied for this interaction or not
    if (mlinfo->binning == true){
        //get new total number of beads and pairs in calculation
        nbeads_old=mlinfo->MLObject->getBeadNumber();
        nbeads_total=nbeads_old+_nbeads_per_frame;

        //if it is a symmetric threebody interaction, less entries are required for the binning
        if (mlinfo->threebody_symmetric == true){
            //number of triplets is always (mlinfo->num_bins)*(mlinfo->num_theta)*(mlinfo->num_bins+1)/2
            ntriples_total=(mlinfo->num_bins)*(mlinfo->num_theta)*(mlinfo->num_bins+1)/2;
            // if second grid (min_out2 < max_out) add number of triplets of second grid
            if (mlinfo->min_out2 < mlinfo->max_out){
                ntriples_total=(mlinfo->num_bins+mlinfo->num_bins2)*(mlinfo->num_theta)*(mlinfo->num_bins+mlinfo->num_bins2+1)/2;
            }
        }
        //if it is not a symmetric threebody interaction, more entries are required for the binning
        if (mlinfo->threebody_symmetric == false){
            //number of triplets is always (mlinfo->num_bins)*(mlinfo->num_theta)*(mlinfo->num_bins)
            ntriples_total=(mlinfo->num_bins)*(mlinfo->num_theta)*(mlinfo->num_bins);
            // if second grid (min_out2 < max_out) add number of triplets of second grid
            if (mlinfo->min_out2 < mlinfo->max_out){
                ntriples_total=(mlinfo->num_bins+mlinfo->num_bins2)*(mlinfo->num_theta)*(mlinfo->num_bins+mlinfo->num_bins2);
            }
        }


        std::cout << "nbeads_total: " << nbeads_total << ", ntriples_total: " << ntriples_total << std::endl;

        //now resize the ML object
        mlinfo->MLObject->Resize(nbeads_total,ntriples_total*3);

        std::cout << "struct number: " << mlinfo->MLObject->getStructNumber() << std::endl;

        //Evaluate Kernel for first grid
        for (int i=0; i<ntriples_total*9; ++i){
            mlinfo->MLObject->setDescriptor(i,mlinfo->binned_structures(i));
        }
        for (int i=0; i<ntriples_total; ++i){
            mlinfo->MLObject->setDescriptorNumber(i,1);
        }
        for (int i=ntriples_total*9; i<ntriples_total*9*2; ++i){
            mlinfo->MLObject->setDescriptor(i,mlinfo->binned_structures(i-ntriples_total*9));
        }
        for (int i=ntriples_total; i<ntriples_total*2; ++i){
            mlinfo->MLObject->setDescriptorNumber(i,2);
        }
        for (int i=ntriples_total*9*2; i<ntriples_total*9*3; ++i){
            mlinfo->MLObject->setDescriptor(i,mlinfo->binned_structures(i-ntriples_total*9*2));
        }
        for (int i=ntriples_total*2; i<ntriples_total*3; ++i){
            mlinfo->MLObject->setDescriptorNumber(i,3);
        }

        double f_cut1, f_cut2;
        double norm12,norm13,theta,dtheta,rinv,cs,save;
        //double norm23

        //iterate over all beads that are taken into account and count descriptors
        for (int i=0; i<_nbeads_per_frame; ++i){
            for (unsigned int j=0; j<conf->getBead(_bead_to_row(i))->DescriptorsSize(); ++j){
                //store descriptors
                Eigen::Vector3d var1 = conf->getBead(_bead_to_row(i))->getDescriptor1(j);
                Eigen::Vector3d var2 = conf->getBead(_bead_to_row(i))->getDescriptor2(j);
                Eigen::Vector3d var3 = conf->getBead(_bead_to_row(i))->getDescriptor3(j);
                int num = conf->getBead(_bead_to_row(i))->getDescriptorNumber(j);

                norm12 = sqrt(var1.x()*var1.x()+var1.y()*var1.y()+var1.z()*var1.z());
                norm13 = sqrt(var2.x()*var2.x()+var2.y()*var2.y()+var2.z()*var2.z());
                //norm23 = sqrt(var3.x()*var3.x()+var3.y()*var3.y()+var3.z()*var3.z());
                rinv = 1.0/(norm12*norm13);
                cs = (var1.x()*var2.x() + var1.y()*var2.y() + var1.z()*var2.z()) * rinv;
                //compute angle between r12 and r13 in degrees
                theta = acos(cs)*180.0/M_PI;
                dtheta = (180.0-mlinfo->min_theta)/(mlinfo->num_theta);
                
                Eigen::Vector3d vec1_eigen(var1.x(), var1.y(), var1.z());
                Eigen::Vector3d vec2_eigen(var2.x(), var2.y(), var2.z());
                Eigen::Vector3d vec3_eigen(var3.x(), var3.y(), var3.z());

                Eigen::Matrix3d Rxz;
                Eigen::Matrix3d Rz;
                Eigen::Matrix3d Rvx;
                Eigen::Matrix3d Rvz;
                Eigen::Matrix3d Rg;

                //if it is a symmetric threebody interaction, d13 always has to be greater than d12
                if (mlinfo->threebody_symmetric == true){
                    //make sure norm13 is always >= norm12
                    if (norm12 > norm13){
                        save = norm12;
                        norm12 = norm13;
                        norm13 = save;
                    }

                    if (vec1_eigen.norm() <= vec2_eigen.norm()) {
                        Rxz = get_Rxz(vec1_eigen);
                        vec1_eigen = Rxz * vec1_eigen;
                        vec2_eigen = Rxz * vec2_eigen;
                        vec3_eigen = Rxz * vec3_eigen;

                        Rvx = get_Rvx(vec1_eigen);
                        vec1_eigen = Rvx * vec1_eigen;
                        vec2_eigen = Rvx * vec2_eigen;
                        vec3_eigen = Rvx * vec3_eigen;

                        Rz = get_Rz(vec2_eigen);
                        vec1_eigen = Rz * vec1_eigen;
                        vec2_eigen = Rz * vec2_eigen;
                        vec3_eigen = Rz * vec3_eigen;

                        Rvz = get_Rvz(vec2_eigen);
                        vec1_eigen = Rvz * vec1_eigen;
                        vec2_eigen = Rvz * vec2_eigen;
                        vec3_eigen = Rvz * vec3_eigen;

                        Rg =  Rvz * Rz * Rvx * Rxz;
                        //now Rg is the total rotation matrix
                    }
                    if (vec1_eigen.norm() > vec2_eigen.norm()) {
                        Rxz = get_Rxz(vec2_eigen);
                        vec1_eigen = Rxz * vec1_eigen;
                        vec2_eigen = Rxz * vec2_eigen;
                        vec3_eigen = Rxz * vec3_eigen;

                        Rvx = get_Rvx(vec2_eigen);
                        vec1_eigen = Rvx * vec1_eigen;
                        vec2_eigen = Rvx * vec2_eigen;
                        vec3_eigen = Rvx * vec3_eigen;

                        Rz = get_Rz(vec1_eigen);
                        vec1_eigen = Rz * vec1_eigen;
                        vec2_eigen = Rz * vec2_eigen;
                        vec3_eigen = Rz * vec3_eigen;

                        Rvz = get_Rvz(vec1_eigen);
                        vec1_eigen = Rvz * vec1_eigen;
                        vec2_eigen = Rvz * vec2_eigen;
                        vec3_eigen = Rvz * vec3_eigen;

                        Rg =  Rvz * Rz * Rvx * Rxz;
                        //now Rg is the total rotation matrix

                        //exchange num 2<->3 because of "swapping" of r12 and r13
                        int store = num;
                        if (store == 2){
                            num = 3;
                        }
                        if (store == 3){
                            num = 2;
                        }
                    }
                }
                //if it is not a symmetric threebody interaction, d13 can also be smaller than d12
                if (mlinfo->threebody_symmetric == false){
                    Rxz = get_Rxz(vec1_eigen);
                    vec1_eigen = Rxz * vec1_eigen;
                    vec2_eigen = Rxz * vec2_eigen;
                    vec3_eigen = Rxz * vec3_eigen;

                    Rvx = get_Rvx(vec1_eigen);
                    vec1_eigen = Rvx * vec1_eigen;
                    vec2_eigen = Rvx * vec2_eigen;
                    vec3_eigen = Rvx * vec3_eigen;

                    Rz = get_Rz(vec2_eigen);
                    vec1_eigen = Rz * vec1_eigen;
                    vec2_eigen = Rz * vec2_eigen;
                    vec3_eigen = Rz * vec3_eigen;

                    Rvz = get_Rvz(vec2_eigen);
                    vec1_eigen = Rvz * vec1_eigen;
                    vec2_eigen = Rvz * vec2_eigen;
                    vec3_eigen = Rvz * vec3_eigen;

                    Rg =  Rvz * Rz * Rvx * Rxz;
                    //now Rg is the total rotation matrix
                }

                f_cut1 = Calculate_fcut(sqrt(var1.x()*var1.x()+var1.y()*var1.y()+var1.z()*var1.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                f_cut2 = Calculate_fcut(sqrt(var2.x()*var2.x()+var2.y()*var2.y()+var2.z()*var2.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);

                //if it is a symmetric threebody interaction, d13 always has to be greater than d12
                if (mlinfo->threebody_symmetric == true){
                    //set entry in mapping matrix
                    //if binning should be done with applying Gaussian smoothing
                    if(mlinfo->gaussian_smearing == true){
                        double out_12,out_13,theta_grid;
                        int pos=0;
                        double normalize_smearing=0.0;
                        double value=0.0;

                        //if only one grid is used
                        if (mlinfo->min_out2 == mlinfo->max_out){
  		            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<mlinfo->num_bins; ++it){
                                out_13 = out_12;
                                for (int jt=it; jt<mlinfo->num_bins; ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
                                        value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
                                        normalize_smearing+=value;
                                        theta_grid += dtheta;
                                        ++pos;
                                    }
                                    out_13 += mlinfo->dx_bins;
                                }
                                out_12 += mlinfo->dx_bins;
                            }
                            pos=0;
                            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<mlinfo->num_bins; ++it){
                                out_13 = out_12;
                                for (int jt=it; jt<mlinfo->num_bins; ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
                                        value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
	                                if ( normalize_smearing==0.0 ){
		                            value = 0.0;
		                        }
			                if ( normalize_smearing>0.0 ){
                                            value /= normalize_smearing;
		                        }
                                        theta_grid += dtheta;

                                        for (unsigned int k=0; k<3; ++k){
                                            for (unsigned int l=0; l<3; ++l){
                                                if (num == 1){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 2){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 3){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3*2+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                            }
                                        }
                                        ++pos;
                                    }
                                    out_13 += mlinfo->dx_bins;
                                }
                                out_12 += mlinfo->dx_bins;
                            }
		        }

		        //if two grids are used
                        if (mlinfo->min_out2 < mlinfo->max_out){
                            pos=0;
                            normalize_smearing=0.0;
                            value=0.0;

		            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<(mlinfo->num_bins+mlinfo->num_bins2); ++it){
                                out_13 = out_12;
                                for (int jt=it; jt<(mlinfo->num_bins+mlinfo->num_bins2); ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
          	                        if (it<mlinfo->num_bins){
                                            value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
         	                        }
 	                                if (it>=mlinfo->num_bins){
	                                    value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
			                }
                                        if (jt<mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
				        }
                                        if (jt>=mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
				        }
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
                                        normalize_smearing+=value;
                                        theta_grid += dtheta;
                                        ++pos;
                                    }
          	                    if ( jt<(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins;
        	                    }
	                            if ( jt==(mlinfo->num_bins-1) ){
                                        out_13 += 0.5*mlinfo->dx_bins;
                                        out_13 += 0.5*mlinfo->dx_bins2;
                	            }
 	                            if ( jt>(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins2;
	                            }
                                }
          	                if ( it<(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins;
        	                }
	                        if ( it==(mlinfo->num_bins-1) ){
                                    out_12 += 0.5*mlinfo->dx_bins;
                                    out_12 += 0.5*mlinfo->dx_bins2;
                	        }
 	                        if ( it>(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins2;
	                        }
                            }
                            pos=0;
		            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<(mlinfo->num_bins+mlinfo->num_bins2); ++it){
                                out_13 = out_12;
                                for (int jt=it; jt<(mlinfo->num_bins+mlinfo->num_bins2); ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
          	                        if (it<mlinfo->num_bins){
                                            value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
         	                        }
 	                                if (it>=mlinfo->num_bins){
	                                    value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
			                }
                                        if (jt<mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
				        }
                                        if (jt>=mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
				        }
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
	                                if ( normalize_smearing==0.0 ){
		                            value = 0.0;
		                        }
			                if ( normalize_smearing>0.0 ){
                                            value /= normalize_smearing;
		                        }
                                        theta_grid += dtheta;
                                        for (unsigned int k=0; k<3; ++k){
                                            for (unsigned int l=0; l<3; ++l){
                                                if (num == 1){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 2){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 3){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3*2+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                            }
                                        }
                                        ++pos;
                                    }
          	                    if ( jt<(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins;
        	                    }
	                            if ( jt==(mlinfo->num_bins-1) ){
                                        out_13 += 0.5*mlinfo->dx_bins;
                                        out_13 += 0.5*mlinfo->dx_bins2;
                	            }
 	                            if ( jt>(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins2;
	                            }
                                }
          	                if ( it<(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins;
        	                }
	                        if ( it==(mlinfo->num_bins-1) ){
                                    out_12 += 0.5*mlinfo->dx_bins;
                                    out_12 += 0.5*mlinfo->dx_bins2;
                	        }
 	                        if ( it>(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins2;
	                        }
                            }
		        }
                    }
                    //if binning should be done without applying Gaussian smoothing
                    if(mlinfo->gaussian_smearing == false){
                        int nr12,nr13,ntheta;
		        //initialization
		        nr12 = 0;
		        nr13 = 0;
		        ntheta = 0;
                        int position;
                        position=0;
                        //if only one grid is used
                        if (mlinfo->min_out2 == mlinfo->max_out){
                            nr12 = (norm12 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                            if (norm12 == (mlinfo->min_out)){
                                nr12 = 0;
                            }
                            nr13 = (norm13 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                            if (norm13 == (mlinfo->min_out)){
                                nr13 = 0;
                            }
                            nr13 -= nr12;
                            ntheta = (theta-mlinfo->min_theta-0.00000001)/dtheta;            
                            if (theta == 180.0){
                                ntheta = (mlinfo->num_bins*2)-1;
                            }
		            if (ntheta < 0){
                                ntheta = 0;
                            }
                            position = 0;
                            for (int in=0; in<nr12; in++){
                                position += (mlinfo->num_bins-in);
                            }
                            position += nr13;                        
                            position *= (mlinfo->num_bins*2);
                            position += ntheta;                 
                        }
                        //if two grids are used, calculate also position2 and position2_1
                        if (mlinfo->min_out2 < mlinfo->max_out){

		            if (norm12 < mlinfo->min_out2){
	                        nr12 = (norm12 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                                if (norm12 == (mlinfo->min_out)){
                                    nr12 = 0;
                                }
                            }
		            if (norm12 >= mlinfo->min_out2){
		                nr12 = mlinfo->num_bins+(norm12 - mlinfo->min_out2 - 0.00000001)/mlinfo->dx_bins2;
                                if (norm12 == (mlinfo->min_out2)){
                                    nr12 = mlinfo->num_bins;
                                }
                            }
		            if (norm13 < mlinfo->min_out2){
	                        nr13 = (norm13 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                                if (norm13 == (mlinfo->min_out)){
                                    nr13 = 0;
                                }
                            }
		            if (norm13 >= mlinfo->min_out2){
		                nr13 = mlinfo->num_bins+(norm13 - mlinfo->min_out2 - 0.00000001)/mlinfo->dx_bins2;
                                if (norm13 == (mlinfo->min_out2)){
                                    nr13 = mlinfo->num_bins;
                                }
                            }
                            nr13 -= nr12;
                            ntheta = (theta-mlinfo->min_theta-0.00000001)/dtheta;
                            if (theta == 180.0){
                                ntheta = ((mlinfo->num_bins+mlinfo->num_bins2)*2)-1;
                            }
		            if (ntheta < 0){
                                ntheta = 0;
                            }
                            position = 0;
                            for (int in=0; in<nr12; in++){
                                position += ((mlinfo->num_bins+mlinfo->num_bins2)-in);
                            }
                            position += nr13;
                            position *= ((mlinfo->num_bins+mlinfo->num_bins2)*2);
                            position += ntheta;
                        }
                        for (unsigned int k=0; k<3; ++k){
                            for (unsigned int l=0; l<3; ++l){
                                if (num == 1){
                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,position*3+l,f_cut1*f_cut2*Rg(l,k));
                                }
                                if (num == 2){
                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,position*3+ntriples_total*3+l,f_cut1*f_cut2*Rg(l,k));
                                }
                                if (num == 3){
                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,position*3+ntriples_total*3*2+l,f_cut1*f_cut2*Rg(l,k));
                                }
                            }
                        }
                    }
                }
                //if it is not a symmetric threebody interaction, d13 can also be smaller than d12
                if (mlinfo->threebody_symmetric == false){
                    //set entry in mapping matrix
                    //if binning should be done with applying Gaussian smoothing
                    if(mlinfo->gaussian_smearing == true){
                        double out_12,out_13,theta_grid;
                        int pos=0;
                        double normalize_smearing=0.0;
                        double value=0.0;

                        //if only one grid is used
                        if (mlinfo->min_out2 == mlinfo->max_out){
  		            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<mlinfo->num_bins; ++it){
                                out_13 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                                for (int jt=0; jt<mlinfo->num_bins; ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
                                        value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
                                        normalize_smearing+=value;
                                        theta_grid += dtheta;
                                        ++pos;
                                    }
                                    out_13 += mlinfo->dx_bins;
                                }
                                out_12 += mlinfo->dx_bins;
                            }
                            pos=0;
                            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<mlinfo->num_bins; ++it){
                                out_13 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                                for (int jt=0; jt<mlinfo->num_bins; ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
                                        value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
	                                if ( normalize_smearing==0.0 ){
		                            value = 0.0;
		                        }
			                if ( normalize_smearing>0.0 ){
                                            value /= normalize_smearing;
		                        }
                                        theta_grid += dtheta;

                                        for (unsigned int k=0; k<3; ++k){
                                            for (unsigned int l=0; l<3; ++l){
                                                if (num == 1){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 2){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 3){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3*2+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                            }
                                        }
                                        ++pos;
                                    }
                                    out_13 += mlinfo->dx_bins;
                                }
                                out_12 += mlinfo->dx_bins;
                            }
		        }

		        //if two grids are used
                        if (mlinfo->min_out2 < mlinfo->max_out){
                            pos=0;
                            normalize_smearing=0.0;
                            value=0.0;

		            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<(mlinfo->num_bins+mlinfo->num_bins2); ++it){
                                out_13 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                                for (int jt=0; jt<(mlinfo->num_bins+mlinfo->num_bins2); ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
          	                        if (it<mlinfo->num_bins){
                                            value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
         	                        }
 	                                if (it>=mlinfo->num_bins){
	                                    value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
			                }
                                        if (jt<mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
				        }
                                        if (jt>=mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
				        }
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
                                        normalize_smearing+=value;
                                        theta_grid += dtheta;
                                        ++pos;
                                    }
          	                    if ( jt<(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins;
        	                    }
	                            if ( jt==(mlinfo->num_bins-1) ){
                                        out_13 += 0.5*mlinfo->dx_bins;
                                        out_13 += 0.5*mlinfo->dx_bins2;
                	            }
 	                            if ( jt>(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins2;
	                            }
                                }
          	                if ( it<(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins;
        	                }
	                        if ( it==(mlinfo->num_bins-1) ){
                                    out_12 += 0.5*mlinfo->dx_bins;
                                    out_12 += 0.5*mlinfo->dx_bins2;
                	        }
 	                        if ( it>(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins2;
	                        }
                            }
                            pos=0;
		            out_12 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                            for (int it=0; it<(mlinfo->num_bins+mlinfo->num_bins2); ++it){
                                out_13 = mlinfo->min_out+0.5*mlinfo->dx_bins;
                                for (int jt=0; jt<(mlinfo->num_bins+mlinfo->num_bins2); ++jt){
                                    theta_grid = mlinfo->min_theta+0.5*dtheta;
                                    for (int kt=0; kt<mlinfo->num_theta; ++kt){
          	                        if (it<mlinfo->num_bins){
                                            value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
         	                        }
 	                                if (it>=mlinfo->num_bins){
	                                    value=exp( - ( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_12 - norm12)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
			                }
                                        if (jt<mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins/mlinfo->smear_scale) ) );
				        }
                                        if (jt>=mlinfo->num_bins){
                                            value*=exp( - ( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) )*( (out_13 - norm13)/(mlinfo->dx_bins2/mlinfo->smear_scale) ) );
				        }
                                        value*=exp( - ( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) )*( (theta_grid-theta)/(dtheta/mlinfo->smear_scale) ) );
                                        //if value is below threshold of 1e-3, replace it with 0.0
                                        if ( value<=0.001 ){
                                            value = 0.0;
                                        }
	                                if ( normalize_smearing==0.0 ){
		                            value = 0.0;
		                        }
			                if ( normalize_smearing>0.0 ){
                                            value /= normalize_smearing;
		                        }
                                        theta_grid += dtheta;
                                        for (unsigned int k=0; k<3; ++k){
                                            for (unsigned int l=0; l<3; ++l){
                                                if (num == 1){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 2){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                                if (num == 3){
                                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+ntriples_total*3*2+l,value*f_cut1*f_cut2*Rg(l,k));
                                                }
                                            }
                                        }
                                        ++pos;
                                    }
          	                    if ( jt<(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins;
        	                    }
	                            if ( jt==(mlinfo->num_bins-1) ){
                                        out_13 += 0.5*mlinfo->dx_bins;
                                        out_13 += 0.5*mlinfo->dx_bins2;
                	            }
 	                            if ( jt>(mlinfo->num_bins-1) ){
                                        out_13 += mlinfo->dx_bins2;
	                            }
                                }
          	                if ( it<(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins;
        	                }
	                        if ( it==(mlinfo->num_bins-1) ){
                                    out_12 += 0.5*mlinfo->dx_bins;
                                    out_12 += 0.5*mlinfo->dx_bins2;
                	        }
 	                        if ( it>(mlinfo->num_bins-1) ){
                                    out_12 += mlinfo->dx_bins2;
	                        }
                            }
		        }
                    }
                    //if binning should be done without applying Gaussian smoothing
                    if(mlinfo->gaussian_smearing == false){
                        int nr12,nr13,ntheta;
		        //initialization
		        nr12 = 0;
		        nr13 = 0;
		        ntheta = 0;
                        int position;
                        position=0;
                        //if only one grid is used
                        if (mlinfo->min_out2 == mlinfo->max_out){
                            nr12 = (norm12 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                            if (norm12 == (mlinfo->min_out)){
                                nr12 = 0;
                            }
                            nr13 = (norm13 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                            if (norm13 == (mlinfo->min_out)){
                                nr13 = 0;
                            }                            
                            ntheta = (theta-mlinfo->min_theta-0.00000001)/dtheta;            
                            if (theta == 180.0){
                                ntheta = (mlinfo->num_bins*2)-1;
                            }
		            if (ntheta < 0){
                                ntheta = 0;
                            }
                            position = nr12*(mlinfo->num_bins);
                            position += nr13;                        
                            position *= (mlinfo->num_bins*2);
                            position += ntheta;                 
                        }
                        //if two grids are used, calculate also position2 and position2_1
                        if (mlinfo->min_out2 < mlinfo->max_out){

		            if (norm12 < mlinfo->min_out2){
	                        nr12 = (norm12 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                                if (norm12 == (mlinfo->min_out)){
                                    nr12 = 0;
                                }
                            }
		            if (norm12 >= mlinfo->min_out2){
		                nr12 = mlinfo->num_bins+(norm12 - mlinfo->min_out2 - 0.00000001)/mlinfo->dx_bins2;
                                if (norm12 == (mlinfo->min_out2)){
                                    nr12 = mlinfo->num_bins;
                                }
                            }
		            if (norm13 < mlinfo->min_out2){
	                        nr13 = (norm13 - mlinfo->min_out - 0.00000001)/mlinfo->dx_bins;
                                if (norm13 == (mlinfo->min_out)){
                                    nr13 = 0;
                                }
                            }
		            if (norm13 >= mlinfo->min_out2){
		                nr13 = mlinfo->num_bins+(norm13 - mlinfo->min_out2 - 0.00000001)/mlinfo->dx_bins2;
                                if (norm13 == (mlinfo->min_out2)){
                                    nr13 = mlinfo->num_bins;
                                }
                            }
                            ntheta = (theta-mlinfo->min_theta-0.00000001)/dtheta;
                            if (theta == 180.0){
                                ntheta = ((mlinfo->num_bins+mlinfo->num_bins2)*2)-1;
                            }
		            if (ntheta < 0){
                                ntheta = 0;
                            }
                            position = nr12*(mlinfo->num_bins+mlinfo->num_bins2);
                            position += nr13;
                            position *= ((mlinfo->num_bins+mlinfo->num_bins2)*2);
                            position += ntheta;
                        }
                        for (unsigned int k=0; k<3; ++k){
                            for (unsigned int l=0; l<3; ++l){
                                if (num == 1){
                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,position*3+l,f_cut1*f_cut2*Rg(l,k));
                                }
                                if (num == 2){
                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,position*3+ntriples_total*3+l,f_cut1*f_cut2*Rg(l,k));
                                }
                                if (num == 3){
                                    mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,position*3+ntriples_total*3*2+l,f_cut1*f_cut2*Rg(l,k));
                                }
                            }
                        }
                    }
                }

            }
            conf->getBead(_bead_to_row(i))->clearDescriptors();
        }
    }
    //if no binning is applied for this interaction
    else{
        //get new total number of beads and pairs in calculation
        nbeads_old=mlinfo->MLObject->getBeadNumber();
        ntriples_old=mlinfo->MLObject->getStructNumber();
        nbeads_total=nbeads_old+_nbeads_per_frame;
        ntriples_total=ntriples_old;

        //iterate over all beads that are taken into account and count descriptors
        for (int i=0; i<_nbeads_per_frame; ++i){
            ntriples_total+=conf->getBead(_bead_to_row(i))->DescriptorsSize();
        }

        std::cout << "nbeads_total: " << nbeads_total << ", ntriples_total: " << ntriples_total << std::endl;

        //now resize the ML object
        mlinfo->MLObject->Resize(nbeads_total,ntriples_total);

        //counter
        int count=0;

        double f_cut1, f_cut2;

        //iterate over all beads that are taken into account and count descriptors
        for (int i=0; i<_nbeads_per_frame; ++i){
            for (unsigned int j=0; j<conf->getBead(_bead_to_row(i))->DescriptorsSize(); ++j){
                //store descriptors
                Eigen::Vector3d var1 = conf->getBead(_bead_to_row(i))->getDescriptor1(j);
                Eigen::Vector3d var2 = conf->getBead(_bead_to_row(i))->getDescriptor2(j);
                Eigen::Vector3d var3 = conf->getBead(_bead_to_row(i))->getDescriptor3(j);
                int num = conf->getBead(_bead_to_row(i))->getDescriptorNumber(j);
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9,var1.x());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 1,var1.y());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 2,var1.z());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 3,var2.x());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 4,var2.y());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 5,var2.z());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 6,var3.x());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 7,var3.y());
                mlinfo->MLObject->setDescriptor((ntriples_old+count)*9 + 8,var3.z());
                mlinfo->MLObject->setDescriptorNumber(ntriples_old+count,num);
                //set entry in mapping matrix
                f_cut1 = Calculate_fcut(sqrt(var1.x()*var1.x()+var1.y()*var1.y()+var1.z()*var1.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                f_cut2 = Calculate_fcut(sqrt(var2.x()*var2.x()+var2.y()*var2.y()+var2.z()*var2.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                for (int k=0; k<3; ++k){
                    mlinfo->MLObject->setMappingMatrix((nbeads_old+i)*3+k,(ntriples_old+count)*3+k,f_cut1*f_cut2);
                }
                //increment counter
                ++count;
            }
            conf->getBead(_bead_to_row(i))->clearDescriptors();
        }
    }
    delete nb;
}

void CGMachineLearning::EvalNonbondedTrain(Topology *conf, MLInfo *mlinfo)
{
    // generate the neighbour list
    NBList *nb;
    bool gridsearch=false;

    if(_options.exists("cg.nbsearch")) {
        if(_options.get("cg.nbsearch").as<string>() == "grid")
            gridsearch=true;
        else if(_options.get("cg.nbsearch").as<string>() == "simple")
            gridsearch=false;
        else throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
    if(gridsearch)
        nb = new NBListGrid();
    else
        nb = new NBList();
    nb->setCutoff(mlinfo->_options->get("ml.max").as<double>()); // implement different cutoffs for different interactions!

    // generate the bead lists
    BeadList beads1, beads2;
    beads1.Generate(*conf, mlinfo->type1);
    beads2.Generate(*conf, mlinfo->type2);
    // is it same types or different types?
    if (mlinfo->type1 == mlinfo->type2)
        nb->Generate(beads1, true);
    else
        nb->Generate(beads1, beads2, true);

    // iterate over all pairs
    for (BeadPair *pair : *nb) {
        int iatom = pair->first()->getId();
        int jatom = pair->second()->getId();
        Eigen::Vector3d var = pair->r();

        //check if iatom and/or jatom of this pair is list _bead_to_row
        for (int i=0; i<_nbeads_per_frame; ++i){
            if (iatom == _bead_to_row(i)){
                // the second and third descriptor values are dummy values, as well as the number
                conf->getBead(_bead_to_row(i))->addDescriptor(var,var,var,1);
            }
            if (jatom == _bead_to_row(i) ){
                // the second and third descriptor values are dummy values, as well as the number
                conf->getBead(_bead_to_row(i))->addDescriptor((-1.0)*var,(-1.0)*var,(-1.0)*var,1);
            }
        }
    }

    int nbeads_old,npairs_old,nbeads_total,npairs_total;

    //check, if binning is applied for this interaction or not
    if (mlinfo->binning == true){
        //get new total number of beads and pairs in calculation
        nbeads_old=mlinfo->MLObject->getBeadNumber();
        nbeads_total=nbeads_old+_nbeads_per_frame;
        //number of pairs is always number of bins
        npairs_total=mlinfo->num_bins;

        std::cout << "nbeads_total: " << nbeads_total << ", npairs_total: " << npairs_total << std::endl;

        //now resize the ML object
        mlinfo->MLObject->Resize(nbeads_total,npairs_total);

        for (int i=0; i<npairs_total*3; ++i){
            mlinfo->MLObject->setDescriptor(i,mlinfo->binned_structures(i));
        }

        double f_cut;
        double norm;

        //iterate over all beads that are taken into account and count descriptors
        for (int i=0; i<_nbeads_per_frame; ++i){
            for (unsigned int j=0; j<conf->getBead(_bead_to_row(i))->DescriptorsSize(); ++j){
                //store descriptor (pair distance)
                Eigen::Vector3d var = conf->getBead(_bead_to_row(i))->getDescriptor1(j);
                norm = sqrt(var.x()*var.x()+var.y()*var.y()+var.z()*var.z());

                f_cut = Calculate_fcut(norm,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);

                Eigen::Vector3d vec_eigen(var.x(), var.y(), var.z());
                Eigen::Matrix3d Rxz = get_Rxz(vec_eigen);
                Eigen::Vector3d vec_eigen_rotate = Rxz*vec_eigen;
                Eigen::Matrix3d Rvx = get_Rvx(vec_eigen_rotate);
                Eigen::Matrix3d R = Rvx*Rxz;
                vec_eigen_rotate = R*vec_eigen;

                //set entry in mapping matrix
                //if binning should be done with applying Gaussian smoothing
                if(mlinfo->gaussian_smearing == true){
                    double normalize_smearing=0.0;
                    double value=0.0;
                    for (int pos=0; pos<mlinfo->num_bins; ++pos){
                        value=exp( - (( (mlinfo->min_out+mlinfo->dx_bins*(pos+0.5)) - norm )/(mlinfo->dx_bins/mlinfo->smear_scale))*(( (mlinfo->min_out+mlinfo->dx_bins*(pos+0.5)) - norm )/(mlinfo->dx_bins/mlinfo->smear_scale)) );
                        //if value is below threshold of 1e-4, replace it with 0.0
                        if ( value<=0.001 ){
                            value = 0.0;
                        }
                        normalize_smearing+=value;
                    }
                    for (int pos=0; pos<mlinfo->num_bins; ++pos){
                        //evaluate Gaussian
                        value=exp( - (( (mlinfo->min_out+mlinfo->dx_bins*(pos+0.5)) - norm )/(mlinfo->dx_bins/mlinfo->smear_scale))*(( (mlinfo->min_out+mlinfo->dx_bins*(pos+0.5)) - norm )/(mlinfo->dx_bins/mlinfo->smear_scale)) );
                        //if value is below threshold of 1e-4, replace it with 0.0
                        if ( value<=0.001 ){
                            value = 0.0;
                        }
                        value /= normalize_smearing;
                        for (unsigned int k=0; k<3; ++k){
                            for (unsigned int l=0; l<3; ++l){
                                mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,pos*3+l,value*f_cut*R(l,k));
                            }
                        }
                    }
                }

                //if binning should be done without applying Gaussian smoothing
                if(mlinfo->gaussian_smearing == false){
                    int position;
                    position = (norm-mlinfo->min_out)/mlinfo->dx_bins;
                    //prevent memory problem
                    if ( norm == (mlinfo->max_out-mlinfo->min_out) ){
                        position = mlinfo->num_bins-1;
                    }
                    if ( norm <= mlinfo->min_out ){
                        position = 0;
                    }
                    for (unsigned int k=0; k<3; ++k){
                        for (unsigned int l=0; l<3; ++l){
                            mlinfo->MLObject->addMappingMatrix((nbeads_old+i)*3+k,position*3+l,f_cut*R(l,k));
                        }
                    }
                }
            }
            conf->getBead(_bead_to_row(i))->clearDescriptors();
        }
    } else {
        //get new total number of beads and pairs in calculation
        nbeads_old=mlinfo->MLObject->getBeadNumber();
        npairs_old=mlinfo->MLObject->getStructNumber();
        nbeads_total=nbeads_old+_nbeads_per_frame;
        npairs_total=npairs_old;

        //iterate over all beads that are taken into account and count descriptors
        for (int i=0; i<_nbeads_per_frame; ++i){
            npairs_total+=conf->getBead(_bead_to_row(i))->DescriptorsSize();
        }

        std::cout << "nbeads_total: " << nbeads_total << ", npairs_total: " << npairs_total << std::endl;

        //now resize the ML object
        mlinfo->MLObject->Resize(nbeads_total,npairs_total);

        int count=0;
        double f_cut;

        //iterate over all beads that are taken into account and count descriptors
        for (int i=0; i<_nbeads_per_frame; ++i){
            for (unsigned int j=0; j<conf->getBead(_bead_to_row(i))->DescriptorsSize(); ++j){
                //store descriptor (pair distance, stored in vector Descriptor1 of the bead)
                Eigen::Vector3d var = conf->getBead(_bead_to_row(i))->getDescriptor1(j);
                mlinfo->MLObject->setDescriptor((npairs_old+count)*3,var.x());
                mlinfo->MLObject->setDescriptor((npairs_old+count)*3 + 1,var.y());
                mlinfo->MLObject->setDescriptor((npairs_old+count)*3 + 2,var.z());
                //set entry in mapping matrix
                f_cut = Calculate_fcut(sqrt(var.x()*var.x()+var.y()*var.y()+var.z()*var.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                for (int k=0; k<3; ++k){
                    mlinfo->MLObject->setMappingMatrix((nbeads_old+i)*3+k,(npairs_old+count)*3+k,f_cut);
                }
                //increment counter
                ++count;
            }
            conf->getBead(_bead_to_row(i))->clearDescriptors();
        }
    }
    delete nb;
}

void CGMachineLearning::EvalNonbondedTest_Threebody(Topology *conf, MLInfo *mlinfo)
{
    // generate the neighbour list
    NBList_3Body *nb;

    bool gridsearch=false;

    if(_options.exists("cg.nbsearch")) {
        if(_options.get("cg.nbsearch").as<string>() == "grid")
            gridsearch=true;
        else if(_options.get("cg.nbsearch").as<string>() == "simple")
            gridsearch=false;
        else throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
    if(gridsearch)
        nb = new NBListGrid_3Body();
    else
        nb = new NBList_3Body();

    //max is the pair cutoff distance for testing
    nb->setCutoff(mlinfo->_options->get("ml.max").as<double>());

    // generate the bead lists
    BeadList beads1, beads2, beads3;
    beads1.Generate(*conf, mlinfo->type1);
    beads2.Generate(*conf, mlinfo->type2);
    beads3.Generate(*conf, mlinfo->type3);

    // Generate the 3body neighbour lists
    // check if type2 and type3 are the same
    if (mlinfo->type2 == mlinfo->type3) {
        // if then type2 and type1 are the same, all three types are the same
        // use the Generate function for this case
        if (mlinfo->type1 == mlinfo->type2) {
            nb->Generate(beads1, true);
        }
        // else use the Generate function for type2 being equal to type3 (and type1 being different)
        if (mlinfo->type1 != mlinfo->type2) {
            nb->Generate(beads1, beads2, true);
        }
    }
    // If type2 and type3 are not the same, use the Generate function for three different bead types
    // (Even if type1 and type2 or type1 and type3 are the same, the Generate function for two different beadtypes
    // is only applicable for the case that type2 is equal to type3
    if (mlinfo->type2 != mlinfo->type3) {
        nb->Generate(beads1, beads2, beads3, true);
    }

    //iterate over all beads that are taken into account and store descriptors
    for (BeadTriple *triple : *nb) {
        int iatom = triple->bead1()->getId();
        int jatom = triple->bead2()->getId();
        int katom = triple->bead3()->getId();
        Eigen::Vector3d var1 = triple->r12();
        Eigen::Vector3d var2 = triple->r13();
        Eigen::Vector3d var3 = triple->r23();
        double d12 = triple->dist12();
        double d13 = triple->dist13();
        //double d23 = triple->dist23();

        //check if iatom and/or jatom of this pair is list _bead_to_row
        for (int i=0; i<_nbeads_per_frame; ++i){
            if (iatom == _bead_to_row(i)){
                // if type2 is equal to type3 ordering is such that d13 is always greater than d12
                if (mlinfo->threebody_symmetric == true){
                    if (d13 > d12){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,1);
                    }
                    if (d12 > d13){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var2,var1,-var3,1);
                    }
                }
                // if type2 is not equal to type3 d13 can be smaller and greater than d12
                if (mlinfo->threebody_symmetric == false){
                    conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,1);
                }
            }
            if (jatom == _bead_to_row(i)){
                // if type2 is equal to type3 ordering is such that d13 is always greater than d12
                if (mlinfo->threebody_symmetric == true){
                    if (d13 > d12){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,2);
                    }
                    if (d12 > d13){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var2,var1,-var3,3);
                    }
                }
                // if type2 is not equal to type3 d13 can be smaller and greater than d12
                if (mlinfo->threebody_symmetric == false){
                    conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,2);
                }
            }
            if (katom == _bead_to_row(i)){
                // if type2 is equal to type3 ordering is such that d13 is always greater than d12
                if (mlinfo->threebody_symmetric == true){
                    if (d13 > d12){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,3);
                    }
                    if (d12 > d13){
                        conf->getBead(_bead_to_row(i))->addDescriptor(var2,var1,-var3,2);
                    }
                }
                // if type2 is not equal to type3 d13 can be smaller and greater than d12
                if (mlinfo->threebody_symmetric == false){
                    conf->getBead(_bead_to_row(i))->addDescriptor(var1,var2,var3,3);
                }
            }
        }
    }

    int ntriples_total=0;

    //iterate over all beads that are taken into account and count descriptors
    for (int i=0; i<_nbeads_per_frame; ++i)
        ntriples_total+=conf->getBead(_bead_to_row(i))->DescriptorsSize();

    std::cout << "ntriples_total: " << ntriples_total << std::endl;

    /// \brief Vector to store descriptors of the test configuration
    Eigen::VectorXd descriptor_test;
    descriptor_test = Eigen::VectorXd::Zero(ntriples_total*9);
    Eigen::VectorXi descriptor_test_number;
    descriptor_test_number = Eigen::VectorXi::Zero(ntriples_total);    
    /// \brief store the mapping matrix for this specific configuration
    Eigen::MatrixXd mapping_test;
    mapping_test = Eigen::MatrixXd::Zero(_nbeads_per_frame*3,ntriples_total*3);

    //counter
    int count=0;

    double f_cut1, f_cut2;

    //iterate over all beads that are taken into account and count descriptors
    for (int i=0; i<_nbeads_per_frame; ++i){
        for (unsigned int j=0; j<conf->getBead(_bead_to_row(i))->DescriptorsSize(); ++j){
            //store descriptors
            Eigen::Vector3d var1 = conf->getBead(_bead_to_row(i))->getDescriptor1(j);
            Eigen::Vector3d var2 = conf->getBead(_bead_to_row(i))->getDescriptor2(j);
            Eigen::Vector3d var3 = conf->getBead(_bead_to_row(i))->getDescriptor3(j);
            int num = conf->getBead(_bead_to_row(i))->getDescriptorNumber(j);            
            descriptor_test(count*9)=var1.x();
            descriptor_test(count*9 + 1)=var1.y();
            descriptor_test(count*9 + 2)=var1.z();
            descriptor_test(count*9 + 3)=var2.x();
            descriptor_test(count*9 + 4)=var2.y();
            descriptor_test(count*9 + 5)=var2.z();
            descriptor_test(count*9 + 6)=var3.x();
            descriptor_test(count*9 + 7)=var3.y();
            descriptor_test(count*9 + 8)=var3.z();
            descriptor_test_number(count)=num;
            f_cut1 = Calculate_fcut(sqrt(var1.x()*var1.x()+var1.y()*var1.y()+var1.z()*var1.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
            f_cut2 = Calculate_fcut(sqrt(var2.x()*var2.x()+var2.y()*var2.y()+var2.z()*var2.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
            //set entry in mapping matrix
            for (int k=0; k<3; ++k)
                mapping_test(i*3+k,count*3+k)=f_cut1*f_cut2;
            //increment counter
            ++count;
        }
        conf->getBead(_bead_to_row(i))->clearDescriptors();
    }

    /// \brief Vector to store the predictions for the bead forces of the test configuration
    Eigen::VectorXd prediction_beads;
    prediction_beads = Eigen::VectorXd::Zero(_nbeads_per_frame*3);

    //now do the prediction for the bead forces
    mlinfo->MLObject->PredictBead_Threebody(descriptor_test,descriptor_test_number,mapping_test,prediction_beads);

    //initialize _test_error
    Eigen::Vector3d test_error (0,0,0);
    //predict the bead energies and compare them to the read in bead energies
    for (int i=0; i<_nbeads_per_frame; ++i){
    	test_error.x()+=abs( prediction_beads(i*3)-_b(i*3) );
    	test_error.y()+=abs( prediction_beads(i*3 + 1)-_b(i*3 + 1) );
    	test_error.z()+=abs( prediction_beads(i*3 + 2)-_b(i*3 + 2) );
    	std::cout << "i: " << i << ", bead_force_x: " << _b(i*3) << ", predicted bead force in x: " << prediction_beads(i*3) << std::endl;
    	std::cout << "i: " << i << ", bead_force_y: " << _b(i*3 + 1) << ", predicted bead force in y: " << prediction_beads(i*3 + 1) << std::endl;
    	std::cout << "i: " << i << ", bead_force_z: " << _b(i*3 + 2) << ", predicted bead force in z: " << prediction_beads(i*3 + 2) << std::endl;
    }
    //normalize
    test_error/=_nbeads_per_frame;
    std::cout << "Interaction " << mlinfo->MLName << ", average test_error for test frame: " << _frame_counter << ": " << test_error << std::endl;

    //update _test_error_sum
    mlinfo->test_error_sum+=(test_error.x() + test_error.y() + test_error.z())/3;

    delete nb;
}

void CGMachineLearning::EvalNonbondedTest(Topology *conf, MLInfo *mlinfo)
{
    // generate the neighbour list
    NBList *nb;
    bool gridsearch=false;

    if(_options.exists("cg.nbsearch")) {
        if(_options.get("cg.nbsearch").as<string>() == "grid")
            gridsearch=true;
        else if(_options.get("cg.nbsearch").as<string>() == "simple")
            gridsearch=false;
        else throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
    if(gridsearch)
        nb = new NBListGrid();
    else
        nb = new NBList();

    //max is the pair cutoff distance for testing
    nb->setCutoff(mlinfo->_options->get("ml.max").as<double>());

    // generate the bead lists
    BeadList beads1, beads2;
    beads1.Generate(*conf, mlinfo->type1);
    beads2.Generate(*conf, mlinfo->type2);
    // is it same types or different types?
    if (mlinfo->type1 == mlinfo->type2)
        nb->Generate(beads1, true);
    else
        nb->Generate(beads1, beads2, true);

    // iterate over all pairs
    for (BeadPair *pair : *nb) {
        int iatom = pair->first()->getId();
        int jatom = pair->second()->getId();
        Eigen::Vector3d var = pair->r();

        //check if iatom and/or jatom of this pair is list _bead_to_row
        for (int i=0; i<_nbeads_per_frame; ++i){
            if (iatom == _bead_to_row(i)){
                // the second and third descriptor values are dummy values, as well as the number
                conf->getBead(_bead_to_row(i))->addDescriptor(var,var,var,1);
            }
            if (jatom == _bead_to_row(i) ){
                // the second and third descriptor values are dummy values, as well as the number
                conf->getBead(_bead_to_row(i))->addDescriptor((-1.0)*var,(-1.0)*var,(-1.0)*var,1);
            }
        }
    }

    int npairs_total=0;

    //iterate over all beads that are taken into account and count descriptors
    for (int i=0; i<_nbeads_per_frame; ++i)
        npairs_total+=conf->getBead(_bead_to_row(i))->DescriptorsSize();

    std::cout << "npairs_total: " << npairs_total << std::endl;

    /// \brief Vector to store descriptors of the test configuration
    Eigen::VectorXd descriptor_test;
    descriptor_test = Eigen::VectorXd::Zero(npairs_total*3);
    /// \brief store the mapping matrix for this specific configuration
    Eigen::MatrixXd mapping_test;
    mapping_test = Eigen::MatrixXd::Zero(_nbeads_per_frame*3,npairs_total*3);

    //pair counter
    int count=0;
    double f_cut;

    //iterate over all beads that are taken into account and count descriptors
    for (int i=0; i<_nbeads_per_frame; ++i){
        for (unsigned int j=0; j<conf->getBead(_bead_to_row(i))->DescriptorsSize(); ++j){
            //store descriptor (pair distance, stored in vector Descriptor1 of the bead)
            Eigen::Vector3d var = conf->getBead(_bead_to_row(i))->getDescriptor1(j);
            descriptor_test(count*3)=var.x();
            descriptor_test(count*3 + 1)=var.y();
            descriptor_test(count*3 + 2)=var.z();
            //set entry in mapping matrix
            f_cut = Calculate_fcut(sqrt(var.x()*var.x()+var.y()*var.y()+var.z()*var.z()),mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
            for (int k=0; k<3; ++k)
                mapping_test(i*3+k,count*3+k)=f_cut;
            //increment counter
            ++count;
        }
        conf->getBead(_bead_to_row(i))->clearDescriptors();
    }

    /// \brief Vector to store the predictions for the bead forces of the test configuration
    Eigen::VectorXd prediction_beads;
    prediction_beads = Eigen::VectorXd::Zero(_nbeads_per_frame*3);

    //now do the prediction for the bead forces
    mlinfo->MLObject->PredictBead(descriptor_test,mapping_test,prediction_beads);

    //initialize _test_error
    Eigen::Vector3d test_error (0,0,0);
    //predict the bead energies and compare them to the read in bead energies
    for (int i=0; i<_nbeads_per_frame; ++i){
    	test_error.x()+=abs( prediction_beads(i*3)-_b(i*3) );
    	test_error.y()+=abs( prediction_beads(i*3 + 1)-_b(i*3 + 1) );
    	test_error.z()+=abs( prediction_beads(i*3 + 2)-_b(i*3 + 2) );
    	std::cout << "i: " << i << ", bead_force_x: " << _b(i*3) << ", predicted bead force in x: " << prediction_beads(i*3) << std::endl;
    	std::cout << "i: " << i << ", bead_force_y: " << _b(i*3 + 1) << ", predicted bead force in y: " << prediction_beads(i*3 + 1) << std::endl;
    	std::cout << "i: " << i << ", bead_force_z: " << _b(i*3 + 2) << ", predicted bead force in z: " << prediction_beads(i*3 + 2) << std::endl;
    }
    //normalize
    test_error/=_nbeads_per_frame;
    std::cout << "interaction " << mlinfo->MLName << ", average test_error for test frame: " << _frame_counter << ": " << test_error << std::endl;

    //update _test_error_sum
    mlinfo->test_error_sum+=(test_error.x() + test_error.y() + test_error.z())/3;

    delete nb;
}

void CGMachineLearning::AccumulateDataTrain()
{
    //now evaluate the ML model
    std::cout << "Accumulate training data" << std::endl;

    //go through all non-bonded interactions and evaluate L_K_L_T

    //initialize matrices with size of first ml object

    //preliminary: so far ML prediction only works for one interaction at a time!
    _L_K_L_T=Eigen::MatrixXd::Zero(_nframes*_nbeads_per_frame*3,_nframes*_nbeads_per_frame*3);
    _L_K_L_T=_mls.front().MLObject->Evaluate_L_K_L_T();

    //regularize _L_K_L_T
    for (int i=0; i<_nframes*_nbeads_per_frame*3; ++i){
        _L_K_L_T(i,i)+=_lambda;
    }
    //_L_K_L_T is now the regularized  L*_K * L^T with _K being the Kernel matrix

    //Simple Least Squares solving of system of lienar equations
    Eigen::HouseholderQR<Eigen::MatrixXd> dec(_L_K_L_T);
    _x=dec.solve(_b);

    //help vector containing residuals
    Eigen::VectorXd residual=_b-_L_K_L_T*_x;
    int count=0;
    //set the coefficients alpha of the ML Kernels
    for (MLInfo &mlinfo : _mls) {
        for (int i=0; i<_nframes*_nbeads_per_frame*3; ++i){
            //each ML object has exactly _nframes*_nbeads_per_frame coefficients
            mlinfo.MLObject->setCoefficient(i,_x[i+(_nframes*_nbeads_per_frame*count*3)]);
        }
        ++count;
    }
    Eigen::VectorXd xo;
    Eigen::VectorXd result_temporary;
    Eigen::VectorXi xon;
    double out_x;
    double theta;
    double dtheta;
    double f_cut1,f_cut2;

    //make predictions for each Kernel (to be adjusted)
    for (MLInfo &mlinfo : _mls) {
        //check if threebody
        if (mlinfo.threebody) {
            xo = Eigen::VectorXd::Zero(mlinfo.num_out*9);
            result_temporary = Eigen::VectorXd::Zero(mlinfo.num_out*3);
            xon = Eigen::VectorXi::Zero(mlinfo.num_out);
            for (unsigned int i=0; i<5; ++i){
                if (i == 0){
                    theta = 70.0;
                }
                if (i == 1){
                    theta = 90.0;
                }
                if (i == 2){
                    theta = 110.0;
                }
                if (i == 3){
                    theta = 130.0;
                }
                if (i == 4){
                    theta = 150.0;
                }
                for (unsigned int j=0; j<3; ++j){
                    //first screen r12 and r13 keeping the angle constant
                    out_x = mlinfo.min_out;
                    for (int k=0; k<mlinfo.num_out; ++k){
                        xo[k*9] = out_x;
                        xo[k*9 + 1] = 0.0;
                        xo[k*9 + 2] = 0.0;
                        xo[k*9 + 3] = out_x*cos(theta * M_PI / 180.0);
                        xo[k*9 + 4] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 5] = 0.0;
                        xo[k*9 + 6] = out_x*(cos(theta * M_PI / 180.0)-1.0);
                        xo[k*9 + 7] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 8] = 0.0;
                        xon[k] = j+1;
                        out_x += mlinfo.dx_out;
                    }
                    mlinfo.MLObject->PredictStruct_Threebody(xo,xon,result_temporary);
                    out_x = mlinfo.min_out;
                    for (int k=0; k<mlinfo.num_out; ++k){
                        xo[k*9] = out_x;
                        xo[k*9 + 1] = 0.0;
                        xo[k*9 + 2] = 0.0;
                        xo[k*9 + 3] = out_x*cos(theta * M_PI / 180.0);
                        xo[k*9 + 4] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 5] = 0.0;
                        xo[k*9 + 6] = out_x*(cos(theta * M_PI / 180.0)-1.0);
                        xo[k*9 + 7] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 8] = 0.0;
                        f_cut1 = Calculate_fcut(out_x,mlinfo._options->get("ml.max").as<double>(),mlinfo.d);
                        f_cut2 = Calculate_fcut(out_x,mlinfo._options->get("ml.max").as<double>(),mlinfo.d);
                        result_temporary(k*3 + 0) *= f_cut1*f_cut2;
                        result_temporary(k*3 + 1) *= f_cut1*f_cut2;
                        result_temporary(k*3 + 2) *= f_cut1*f_cut2;
                        out_x += mlinfo.dx_out;
                    }
                    for (int k=0; k<mlinfo.num_out; ++k){
                        mlinfo.result(3*i+j,k*3 + 0) = result_temporary(k*3 + 0);
                        mlinfo.result(3*i+j,k*3 + 1) = result_temporary(k*3 + 1);
                        mlinfo.result(3*i+j,k*3 + 2) = result_temporary(k*3 + 2);
                    }                
                }
            }
            dtheta = 180.0/(mlinfo.num_out-1);
            for (unsigned int i=0; i<4; ++i){
                if (i == 0){
                    out_x = mlinfo.min_out;
                }
                if (i == 1){
                    out_x = mlinfo.min_out+(mlinfo.max_out-mlinfo.min_out)*1.0/3.0;
                }
                if (i == 2){
                    out_x = mlinfo.min_out+(mlinfo.max_out-mlinfo.min_out)*2.0/3.0;
                }
                if (i == 3){
                    out_x = mlinfo.max_out;
                }
                for (unsigned int j=0; j<3; ++j){
                    //now screen the angle and keep r12 and r13 constant
                    theta = 0.0;
                    for (int k=0; k<mlinfo.num_out; ++k){
                        xo[k*9] = out_x;
                        xo[k*9 + 1] = 0.0;
                        xo[k*9 + 2] = 0.0;
                        xo[k*9 + 3] = out_x*cos(theta * M_PI / 180.0);
                        xo[k*9 + 4] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 5] = 0.0;
                        xo[k*9 + 6] = out_x*(cos(theta * M_PI / 180.0)-1.0);
                        xo[k*9 + 7] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 8] = 0.0;
                        xon[k] = j+1;
                        theta += dtheta;
                    }
                    mlinfo.MLObject->PredictStruct_Threebody(xo,xon,result_temporary);
                    theta = 0.0;
                    for (int k=0; k<mlinfo.num_out; ++k){
                        xo[k*9] = out_x;
                        xo[k*9 + 1] = 0.0;
                        xo[k*9 + 2] = 0.0;
                        xo[k*9 + 3] = out_x*cos(theta * M_PI / 180.0);
                        xo[k*9 + 4] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 5] = 0.0;
                        xo[k*9 + 6] = out_x*(cos(theta * M_PI / 180.0)-1.0);
                        xo[k*9 + 7] = out_x*sin(theta * M_PI / 180.0);
                        xo[k*9 + 8] = 0.0;
                        f_cut1 = Calculate_fcut(out_x,mlinfo._options->get("ml.max").as<double>(),mlinfo.d);
                        f_cut2 = Calculate_fcut(out_x,mlinfo._options->get("ml.max").as<double>(),mlinfo.d);
                        result_temporary(k*3 + 0) *= f_cut1*f_cut2;
                        result_temporary(k*3 + 1) *= f_cut1*f_cut2;
                        result_temporary(k*3 + 2) *= f_cut1*f_cut2;
                        theta += dtheta;
                    }
                    for (int k=0; k<mlinfo.num_out; ++k){
                        mlinfo.resulttheta(3*i+j,k*3 + 0) = result_temporary(k*3 + 0);
                        mlinfo.resulttheta(3*i+j,k*3 + 1) = result_temporary(k*3 + 1);
                        mlinfo.resulttheta(3*i+j,k*3 + 2) = result_temporary(k*3 + 2);
                    }                
                }
            }
        } else {
            xo = Eigen::VectorXd::Zero(mlinfo.num_out*3);
            result_temporary = Eigen::VectorXd::Zero(mlinfo.num_out*3);
            //screen r
            out_x = mlinfo.min_out;
            for (int i=0; i<mlinfo.num_out; ++i){
                xo[i*3] = out_x; //out_x;
                xo[i*3 + 1] = 0.0; //out_y;
                xo[i*3 + 2] = 0.0; //out_z;
                out_x += mlinfo.dx_out;
            }
            mlinfo.MLObject->PredictStruct(xo,result_temporary);
            out_x = mlinfo.min_out;
            for (int i=0; i<mlinfo.num_out; ++i){
                xo[i*3] = out_x; //out_x;
                xo[i*3 + 1] = 0.0; //out_y;
                xo[i*3 + 2] = 0.0; //out_z;
                f_cut1 = Calculate_fcut(out_x,mlinfo._options->get("ml.max").as<double>(),mlinfo.d);
                result_temporary(i*3 + 0) *= f_cut1;
                result_temporary(i*3 + 1) *= f_cut1;
                result_temporary(i*3 + 2) *= f_cut1;
                out_x += mlinfo.dx_out;
            }
            for (int i=0; i<mlinfo.num_out; ++i){
                mlinfo.result(0,i*3 + 0) = result_temporary(i*3 + 0);
                mlinfo.result(0,i*3 + 1) = result_temporary(i*3 + 1);
                mlinfo.result(0,i*3 + 2) = result_temporary(i*3 + 2);
            }
        }
    }
}

void CGMachineLearning::AccumulateDataTest()
{
    std::cout << "Accumulate test data" << std::endl;

    MLContainer::iterator mliter;
    //normalize test_error_sum for each nonbonded interaction
    for (MLInfo &mlinfo : _mls) {
        mlinfo.test_error_sum/=_frame_counter;
        std::cout << "Average error of interaction " << mlinfo.MLName << " of test set consisting of " << _frame_counter << " frames: " << mlinfo.test_error_sum << std::endl;
    }
}

void CGMachineLearning::EvaluateTable_Threebody(MLInfo *mlinfo)
{
    //output vector, to be changed
    Eigen::VectorXd xo;
    Eigen::VectorXi xon;
    Eigen::VectorXd result_prel;

    double out_12;
    double out_13;
    double theta;
    double dtheta;
    int count;
    double f_cut1,f_cut2;
    
    dtheta = 180.0/(mlinfo->num_table*2);

    // if type2 is equal to type3 ordering is such that d13 is always greater than d12
    if (mlinfo->threebody_symmetric == true){
        //make predictions for each Kernel (to be adjusted)
        xo = Eigen::VectorXd::Zero(mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1)*9);
        result_prel = Eigen::VectorXd::Zero(mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1)*3);
        xon = Eigen::VectorXi::Zero(mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1));
        std::cout << "Evaluating Table for threebody interaction: " << mlinfo->MLName << std::endl;

        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = out_12;
            for (int j=i; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    xon[count] = 1;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        mlinfo->MLObject->PredictStruct_Threebody(xo,xon,result_prel);
        //now multiply the predicted forces with the cutoff function
        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = out_12;
            for (int j=i; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    f_cut1 = Calculate_fcut(out_12,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    f_cut2 = Calculate_fcut(out_13,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    result_prel[count*3 + 0] *= f_cut1*f_cut2;
                    result_prel[count*3 + 1] *= f_cut1*f_cut2;
                    result_prel[count*3 + 2] *= f_cut1*f_cut2;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        for (int i=0; i<( mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1) ); ++i){
            mlinfo->output1_1[i] = ( (result_prel[i*3 + 0]*xo[i*9 + 4] - result_prel[i*3 + 1]*xo[i*9 + 3])/(xo[i*9 + 4]*xo[i*9 + 0] - xo[i*9 + 3]*xo[i*9 + 1]) );
            mlinfo->output1_2[i] = ( (result_prel[i*3 + 1]*xo[i*9 + 0] - result_prel[i*3 + 0]*xo[i*9 + 1])/(xo[i*9 + 4]*xo[i*9 + 0] - xo[i*9 + 3]*xo[i*9 + 1]) );
        }

        for (int i=0; i<( mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1) ); ++i){
            xon[i] = 2;
        }
        mlinfo->MLObject->PredictStruct_Threebody(xo,xon,result_prel);
        //now multiply the predicted forces with the cutoff function
        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = out_12;
            for (int j=i; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    f_cut1 = Calculate_fcut(out_12,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    f_cut2 = Calculate_fcut(out_13,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    result_prel[count*3 + 0] *= f_cut1*f_cut2;
                    result_prel[count*3 + 1] *= f_cut1*f_cut2;
                    result_prel[count*3 + 2] *= f_cut1*f_cut2;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        for (int i=0; i<( mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1) ); ++i){
            mlinfo->output2_1[i] = ( (result_prel[i*3 + 0]*xo[i*9 + 7] - result_prel[i*3 + 1]*xo[i*9 + 6])/(xo[i*9 + 7]*xo[i*9 + 0] - xo[i*9 + 6]*xo[i*9 + 1]) );
            mlinfo->output2_2[i] = ( (result_prel[i*3 + 1]*xo[i*9 + 0] - result_prel[i*3 + 0]*xo[i*9 + 1])/(xo[i*9 + 7]*xo[i*9 + 0] - xo[i*9 + 6]*xo[i*9 + 1]) );
        }

        for (int i=0; i<( mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1) ); ++i){
            xon[i] = 3;
        }
        mlinfo->MLObject->PredictStruct_Threebody(xo,xon,result_prel);
        //now multiply the predicted forces with the cutoff function
        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = out_12;
            for (int j=i; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    f_cut1 = Calculate_fcut(out_12,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    f_cut2 = Calculate_fcut(out_13,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    result_prel[count*3 + 0] *= f_cut1*f_cut2;
                    result_prel[count*3 + 1] *= f_cut1*f_cut2;
                    result_prel[count*3 + 2] *= f_cut1*f_cut2;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        for (int i=0; i<( mlinfo->num_table*mlinfo->num_table*(mlinfo->num_table+1) ); ++i){
            mlinfo->output3_1[i] = ( (result_prel[i*3 + 0]*xo[i*9 + 7] - result_prel[i*3 + 1]*xo[i*9 + 6])/(xo[i*9 + 7]*xo[i*9 + 3] - xo[i*9 + 6]*xo[i*9 + 4]) );
            mlinfo->output3_2[i] = ( (result_prel[i*3 + 1]*xo[i*9 + 3] - result_prel[i*3 + 0]*xo[i*9 + 4])/(xo[i*9 + 7]*xo[i*9 + 3] - xo[i*9 + 6]*xo[i*9 + 4]) );
        }
    }
    // if type2 is not equal to type3 ordering is such that d13 can also be smaller than d12
    if (mlinfo->threebody_symmetric == false){
        //make predictions for each Kernel (to be adjusted)
        xo = Eigen::VectorXd::Zero(2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table*9);
        result_prel = Eigen::VectorXd::Zero(2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table*3);
        xon = Eigen::VectorXi::Zero(2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table);
        std::cout << "Evaluating Table for threebody interaction: " << mlinfo->MLName << std::endl;

        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = mlinfo->min_out+0.5*mlinfo->dx_table;
            for (int j=0; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    xon[count] = 1;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        mlinfo->MLObject->PredictStruct_Threebody(xo,xon,result_prel);
        //now multiply the predicted forces with the cutoff function
        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = mlinfo->min_out+0.5*mlinfo->dx_table;
            for (int j=0; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    f_cut1 = Calculate_fcut(out_12,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    f_cut2 = Calculate_fcut(out_13,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    result_prel[count*3 + 0] *= f_cut1*f_cut2;
                    result_prel[count*3 + 1] *= f_cut1*f_cut2;
                    result_prel[count*3 + 2] *= f_cut1*f_cut2;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        for (int i=0; i<( 2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table ); ++i){
            mlinfo->output1_1[i] = ( (result_prel[i*3 + 0]*xo[i*9 + 4] - result_prel[i*3 + 1]*xo[i*9 + 3])/(xo[i*9 + 4]*xo[i*9 + 0] - xo[i*9 + 3]*xo[i*9 + 1]) );
            mlinfo->output1_2[i] = ( (result_prel[i*3 + 1]*xo[i*9 + 0] - result_prel[i*3 + 0]*xo[i*9 + 1])/(xo[i*9 + 4]*xo[i*9 + 0] - xo[i*9 + 3]*xo[i*9 + 1]) );
        }

        for (int i=0; i<( 2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table ); ++i){
            xon[i] = 2;
        }
        mlinfo->MLObject->PredictStruct_Threebody(xo,xon,result_prel);
        //now multiply the predicted forces with the cutoff function
        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = mlinfo->min_out+0.5*mlinfo->dx_table;
            for (int j=0; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    f_cut1 = Calculate_fcut(out_12,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    f_cut2 = Calculate_fcut(out_13,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    result_prel[count*3 + 0] *= f_cut1*f_cut2;
                    result_prel[count*3 + 1] *= f_cut1*f_cut2;
                    result_prel[count*3 + 2] *= f_cut1*f_cut2;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        for (int i=0; i<( 2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table ); ++i){
            mlinfo->output2_1[i] = ( (result_prel[i*3 + 0]*xo[i*9 + 7] - result_prel[i*3 + 1]*xo[i*9 + 6])/(xo[i*9 + 7]*xo[i*9 + 0] - xo[i*9 + 6]*xo[i*9 + 1]) );
            mlinfo->output2_2[i] = ( (result_prel[i*3 + 1]*xo[i*9 + 0] - result_prel[i*3 + 0]*xo[i*9 + 1])/(xo[i*9 + 7]*xo[i*9 + 0] - xo[i*9 + 6]*xo[i*9 + 1]) );
        }

        for (int i=0; i<( 2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table ); ++i){
            xon[i] = 3;
        }
        mlinfo->MLObject->PredictStruct_Threebody(xo,xon,result_prel);
        //now multiply the predicted forces with the cutoff function
        out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
        theta = 0.0+0.5*dtheta;
        count = 0;
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = mlinfo->min_out+0.5*mlinfo->dx_table;
            for (int j=0; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    xo[count*9 + 0] = out_12;
                    xo[count*9 + 1] = 0.0;
                    xo[count*9 + 2] = 0.0;
                    xo[count*9 + 3] = out_13*cos(theta * M_PI / 180.0);
                    xo[count*9 + 4] = out_13*sin(theta * M_PI / 180.0);
                    xo[count*9 + 5] = 0.0;
                    xo[count*9 + 6] = xo[count*9 + 3] - xo[count*9 + 0];
                    xo[count*9 + 7] = xo[count*9 + 4]; //y component of r13 and r23 are the same!!
                    xo[count*9 + 8] = 0.0;
                    f_cut1 = Calculate_fcut(out_12,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    f_cut2 = Calculate_fcut(out_13,mlinfo->_options->get("ml.max").as<double>(),mlinfo->d);
                    result_prel[count*3 + 0] *= f_cut1*f_cut2;
                    result_prel[count*3 + 1] *= f_cut1*f_cut2;
                    result_prel[count*3 + 2] *= f_cut1*f_cut2;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }

        for (int i=0; i<( 2*mlinfo->num_table*mlinfo->num_table*mlinfo->num_table ); ++i){
            mlinfo->output3_1[i] = ( (result_prel[i*3 + 0]*xo[i*9 + 7] - result_prel[i*3 + 1]*xo[i*9 + 6])/(xo[i*9 + 7]*xo[i*9 + 3] - xo[i*9 + 6]*xo[i*9 + 4]) );
            mlinfo->output3_2[i] = ( (result_prel[i*3 + 1]*xo[i*9 + 3] - result_prel[i*3 + 0]*xo[i*9 + 4])/(xo[i*9 + 7]*xo[i*9 + 3] - xo[i*9 + 6]*xo[i*9 + 4]) );
        }
    }
}

void CGMachineLearning::WriteTable_Threebody(MLInfo *mlinfo)
{
    std::cout << "Writing out table for threebody interaction: " << mlinfo->MLName << std::endl;

    string file_extension;
    string file_name;

    double out_12;
    double out_13;
    double theta;
    double dtheta;
    int count;

    ofstream out;

    // construct meaningful outfile name
    file_name = mlinfo->MLName;
    file_extension = ".table";        
    file_name = file_name + file_extension;

    // print output file names on stdout
    std::cout << "Updating file: " << file_name << std::endl;

    out.open(file_name.c_str());
    out << "ENTRY1" << std::endl;
    out << "N " << mlinfo->num_table << " rmin " << mlinfo->min_out+0.5*mlinfo->dx_table << " rmax " << mlinfo->max_out-0.5*mlinfo->dx_table << std::endl;
    out << std::endl;

    out_12 = mlinfo->min_out+0.5*mlinfo->dx_table;
    dtheta = 180.0/(mlinfo->num_table*2);
    theta = 0.0+0.5*dtheta;
    count = 0;

    // if type2 is equal to type3 ordering is such that d13 is always greater than d12
    if (mlinfo->threebody_symmetric == true){
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = out_12;
            for (int j=i; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    out << count+1 << " " << out_12 << " " << out_13 << " " << theta << " " << mlinfo->output1_1[count] << " " << mlinfo->output1_2[count] << " " << mlinfo->output2_1[count] << " " << mlinfo->output2_2[count] << " " << mlinfo->output3_1[count] << " " << mlinfo->output3_2[count] << " " << " 0.0" << endl;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }
    }
    // if type2 is not equal to type3 ordering is such that d13 can also be smaller than d12
    if (mlinfo->threebody_symmetric == false){
        for (int i=0; i<mlinfo->num_table; ++i){
            out_13 = mlinfo->min_out+0.5*mlinfo->dx_table;
            for (int j=0; j<mlinfo->num_table; ++j){
                theta = 0.0+0.5*dtheta;
                for (int k=0; k<(mlinfo->num_table*2); ++k){
                    out << count+1 << " " << out_12 << " " << out_13 << " " << theta << " " << mlinfo->output1_1[count] << " " << mlinfo->output1_2[count] << " " << mlinfo->output2_1[count] << " " << mlinfo->output2_2[count] << " " << mlinfo->output3_1[count] << " " << mlinfo->output3_2[count] << " " << " 0.0" << endl;
                    theta += dtheta;
                    ++count;
                }
                out_13 += mlinfo->dx_table;
            }
            out_12 += mlinfo->dx_table;
        }
    }

    out.close();
}

const double CGMachineLearning::Calculate_fcut(double r, double r_cut, double d)
{
    double result;
    if (r <= (r_cut-d)){
        result = 1.0;
    }
    if ( (r > (r_cut-d)) && (r < (r_cut)) ){
        result = ( cos( M_PI*( (r-r_cut+d)/d ) ) + 1 )/2.0;
    }
    if ( r >= r_cut ){
        result = 0.0;
    }
    return result;
}

const Eigen::Matrix3d CGMachineLearning::get_Rxz(Eigen::Vector3d &vec)
{
    double theta = atan(vec(0) / vec(1));
    double alpha = atan( ( vec(0) * sin(theta) + vec(1) * cos(theta) ) / vec(2) );

    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();

    double ct = cos(theta);
    double st = sin(theta);
    double ca = cos(alpha);
    double sa = sin(alpha);

    R(0,0) = ct;
    R(1,0) = ca * st;
    R(2,0) = sa * st;
    R(0,1) = -st;
    R(1,1) = ca * ct;
    R(2,1) = sa * ct;
    R(1,2) = -sa;
    R(2,2) = ca;

    return R;
}

const Eigen::Matrix3d CGMachineLearning::get_Rz(Eigen::Vector3d &vec)
{
    double theta = atan(vec(0) / vec(1));

    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();

    double ct = cos(theta);
    double st = sin(theta);

    R(0,0) = ct;
    R(1,0) = st;
    R(0,1) = -st;
    R(1,1) = ct;
    R(2,2) = 1;

    return R;
}

const Eigen::Matrix3d CGMachineLearning::get_Rvx(Eigen::Vector3d &vec)
{
    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
    R(0,0) = 1;
    R(1,1) = 1;
    R(2,2) = 1;

    if (vec(2) < 0) {
    	R(1,1) = -1;
    	R(2,2) = -1;
    }
    return R;
}

const Eigen::Matrix3d CGMachineLearning::get_Rvz(Eigen::Vector3d &vec)
{
    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
    R(0,0) = 1;
    R(1,1) = 1;
    R(2,2) = 1;

    if (vec(1) < 0) {
    	R(0,0) = -1;
    	R(1,1) = -1;
    }
    return R;
}

void CGMachineLearning::LoadOptions(const string &file) 
{
  _options.LoadFromXML(file);
  _nonbonded = _options.Select("cg.non-bonded");
}
