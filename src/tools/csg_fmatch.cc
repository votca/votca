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

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <votca/tools/cubicspline.h>
#include <votca/csg/nblistgrid.h>
#include <votca/csg/beadlist.h>
#include "csg_fmatch.h"
#include <votca/tools/table.h>
#include <votca/tools/linalg.h>

int main(int argc, char** argv)
{
    CGForceMatching app;
    return app.Exec(argc, argv);
}


void CGForceMatching::Initialize(void)
{
    CsgApplication::Initialize();
    AddProgramOptions()
        ("options", boost::program_options::value<string>(), "  options file for coarse graining")
        ("trj-force", boost::program_options::value<string>(), "  coarse-grained trajectory containing forces of already known interactions");
}

bool CGForceMatching::EvaluateOptions()
{
    CsgApplication::EvaluateOptions();
    CheckRequired("trj", "no trajectory file specified");
    CheckRequired("options", "need to specify options file");
    LoadOptions(OptionsMap()["options"].as<string>());

    _has_existing_forces = false;
    if(OptionsMap().count("trj-force"))
        _has_existing_forces = true;
    return true;
}

void CGForceMatching::BeginEvaluate(Topology *top, Topology *top_atom)
{
    // set counters to zero value:
    _nblocks = 0;
    _line_cntr = _col_cntr = 0;

    // Number of CG beads in topology
    _nbeads = top->BeadCount();
    // Set frame counter to zero
    _frame_counter = 0;

    // accuracy for evaluating the difference in bead positions (default 1e-5)
    _dist = 1e-5;
    if (_options.exists("cg.fmatch.dist")) {
        _dist = _options.get("cg.fmatch.dist").as<double>();
    }
    //std::cout << "_dist: " << _dist << std::endl;  
    
    // read _nframes from input file
    _nframes = _options.get("cg.fmatch.frames_per_block").as<int>();
    // read _constr_least_sq from input file
    _constr_least_sq = _options.get("cg.fmatch.constrainedLS").as<bool>();
        
    // initializing bonded interactions
    for (list<Property*>::iterator iter = _bonded.begin();
            iter != _bonded.end(); ++iter) {
        SplineInfo *i = new SplineInfo(_splines.size(), true, _col_cntr, *iter);
        //adjust initial matrix dimensions:
        _line_cntr += i->num_gridpoints;
        _col_cntr += 2 * i->num_gridpoints;
        //if periodic potential, one additional constraint has to be taken into account -> 1 additional line in matrix
        if(i->periodic != 0){
            _line_cntr += 1;
        }
       
        // add spline to container
        _splines.push_back(i);
    }

    // initializing non-bonded interactions
    for (list<Property*>::iterator iter = _nonbonded.begin();
            iter != _nonbonded.end(); ++iter) {
        SplineInfo *i = new SplineInfo(_splines.size(), false, _col_cntr, *iter);
        //adjust initial matrix dimensions:
        _line_cntr += i->num_gridpoints;
        _col_cntr += 2 * i->num_gridpoints;

        // add spline to container
        _splines.push_back(i);
    }

    cout << "\nYou are using VOTCA!\n";
    cout << "\nhey, somebody wants to forcematch!\n";


    // now initialize _A, _b, _x and probably _B_constr
    // depending on least-squares algorithm used
    if (_constr_least_sq) { // Constrained Least Squares
        
        cout << "\nUsing constrained Least Squares!\n " << endl;

        // assign _least_sq_offset
        _least_sq_offset = 0;

        // resize and clear _B_constr
        _B_constr=Eigen::MatrixXd::Zero(_line_cntr, _col_cntr);

        // resize matrix _A
        _A=Eigen::MatrixXd::Zero(3 * _nbeads *_nframes, _col_cntr);
        // resize vector _b
        _b=Eigen::VectorXd::Zero(3 * _nbeads *_nframes);

        // in case of constrained least squares smoothing conditions
        // are assigned to matrix _B_constr
        FmatchAssignSmoothCondsToMatrix(_B_constr);
    } else { // Simple Least Squares

        cout << "\nUsing simple Least Squares! " << endl;
        // assign _least_sq_offset
        _least_sq_offset = _line_cntr;

        // resize matrix _A
        _A=Eigen::MatrixXd::Zero(_line_cntr + 3 * _nbeads *_nframes, _col_cntr);
        // resize vector _b
        _b=Eigen::VectorXd::Zero(_line_cntr + 3 * _nbeads *_nframes); 

        // in case of simple least squares smoothing conditions
        // are assigned to matrix _A
        FmatchAssignSmoothCondsToMatrix(_A);
        // clear _b (only necessary in simple least squares)
        _b.setZero();
    }
    // resize and clear _x
    _x=Eigen::VectorXd::Zero(_col_cntr);

    if(_has_existing_forces) {
        _top_force.CopyTopologyData(top);
        _trjreader_force = TrjReaderFactory().Create(_op_vm["trj-force"].as<string>());
        if(_trjreader_force == NULL)
            throw runtime_error(string("input format not supported: ") + _op_vm["trj-force"].as<string>());
        // open the trajectory
        _trjreader_force->Open(_op_vm["trj-force"].as<string>());
        // read in first frame
        _trjreader_force->FirstFrame(_top_force);
    }
}

CGForceMatching::SplineInfo::SplineInfo(int index, bool bonded_, int matr_pos_, Property *options) 
{
    // initialize standard data
    splineIndex = index;
    _options = options;
    splineName = options->get("name").value();
    bonded = bonded_;
    //in general natural boundary conditions are used for splines
    periodic = 0;

    // get non-bonded information
    if (!bonded) {
        type1 = options->get("type1").value();
        type2 = options->get("type2").value(); 
    }
    if (bonded) {
        //check if option periodic exists
        if (options->exists("fmatch.periodic")) {
            periodic = options->get("fmatch.periodic").as<bool>();
            // set cubic spline Spline boundary conditions to periodic 
            Spline.setBCInt(1);            
        }        
    }

    // initialize the grid
    double grid_min = options->get("fmatch.min").as<double>();
    double grid_max = options->get("fmatch.max").as<double>();
    double grid_step = options->get("fmatch.step").as<double>();

    // GenerateGrid returns number of grid points. We subtract 1 to get
    // the number of spline functions
    num_gridpoints = Spline.GenerateGrid(grid_min, grid_max, grid_step);
    num_splinefun = num_gridpoints - 1;

    cout << "Number of spline functions for the interaction " << splineName << ":" << num_splinefun << endl;

    matr_pos = matr_pos_;

    // initialize grid for block averaging
    dx_out = options->get("fmatch.out_step").as<double>();
    // number of output grid points
    num_outgrid = 1 + (int)((grid_max-grid_min)/dx_out);
    result=Eigen::VectorXd::Zero(num_outgrid);
    error=Eigen::VectorXd::Zero(num_outgrid);
    resSum=Eigen::VectorXd::Zero(num_outgrid);
    resSum2=Eigen::VectorXd::Zero(num_outgrid);
    block_res_f=Eigen::VectorXd::Zero(num_outgrid);
    block_res_f2=Eigen::VectorXd::Zero(num_outgrid);

}
void CGForceMatching::EndEvaluate()
{
    // sanity check
    if (_nblocks == 0) {
        cerr << "\nERROR in CGForceMatching::EndEvaluate - No blocks have been processed so far" << endl;
        cerr << "It might be that you are using trajectory, which is smaller than needed for one block" << endl;
        cerr << "Check your input!" << endl;
        exit(-1);
    }
     
    cout << "\nWe are done, thank you very much!" << endl;
    if(_has_existing_forces) {
        _trjreader_force->Close();
        delete _trjreader_force;
    }
}

void CGForceMatching::WriteOutFiles()
{
    string file_extension = ".force";
    string file_name;
    Table force_tab;

    // table with error column
    force_tab.SetHasYErr(true);

    SplineContainer::iterator is;

    for (is = _splines.begin(); is != _splines.end(); ++is) {
        // construct meaningful outfile name
        file_name = (*is)->splineName;
        file_name = file_name + file_extension;
        
        // resize table
        force_tab.resize((*is)->num_outgrid);

        // print output file names on stdout
        cout << "Updating file: " << file_name << endl;

        // loop over output grid points
        for (int i = 0; i < (*is)->num_outgrid; i++) {
            // average value
            (*is)->result[i] = (*is)->resSum[i] / _nblocks;
            // standard deviation of the average
            (*is)->error[i] = sqrt( ((*is)->resSum2[i] / _nblocks - (*is)->result[i] * (*is)->result[i])/_nblocks );
        }

        // first output point = first grid point
        double out_x = (*is)->Spline.getGridPoint(0);
        // loop over output grid
        for (int i = 0; i < (*is)->num_outgrid; i++) {
            // put point, result, flag and accuracy at point out_x into the table
            force_tab.set(i, out_x, (-1.0) * (*is)->result[i], 'i', (*is)->error[i]);
            // update out_x for the next iteration
            out_x += (*is)->dx_out;
        }
        // save table in the file
        force_tab.Save(file_name);
        // clear the table for the next spline
        force_tab.clear();
    }
}

void CGForceMatching::EvalConfiguration(Topology *conf, Topology *conf_atom) 
{
    SplineContainer::iterator spiter;
    if(conf->BeadCount() == 0)
        throw std::runtime_error("CG Topology has 0 beads, check your mapping file!");
    if(_has_existing_forces) {
        if(conf->BeadCount() != _top_force.BeadCount())
            throw std::runtime_error("number of beads in topology and force topology does not match");
        for(int i=0; i<conf->BeadCount(); ++i) {                
//            cout << "conf->getBead(" << i << ")->getPos(): " << conf->getBead(i)->getPos() << endl;                
//            cout << "_top_force->getBead(" << i << ")->getPos(): " << _top_force.getBead(i)->getPos() << endl; 
            conf->getBead(i)->F() -= _top_force.getBead(i)->getF();
            vec d = conf->getBead(i)->getPos() - _top_force.getBead(i)->getPos();
//            cout << "vec d of bead " << i << ": " << d << "abs(d): " << abs(d) << endl; 
            if(abs(d) > _dist)//default is 1e-5, otherwise it can be a too strict criterion
//                cout << "conf->getBead(" << i << ")->getPos(): " << conf->getBead(i)->getPos() << endl;                
//                cout << "_top_force->getBead(" << i << ")->getPos(): " << _top_force.getBead(i)->getPos() << endl; 
                throw std::runtime_error("One or more bead positions in mapped and reference force trajectory differ by more than 1e-5");
        }
    }

    for (spiter = _splines.begin(); spiter != _splines.end(); ++spiter) {
        SplineInfo *sinfo = *spiter;
        if (sinfo->bonded) // bonded interaction
            EvalBonded(conf, sinfo);
        else // non-bonded interaction
            EvalNonbonded(conf, sinfo);
    }
    
    // loop for the forces vector: 
    // hack, chage the Has functions..
    if (conf->getBead(0)->HasF()) {
        vec Force(0., 0., 0.);
        for (int iatom = 0; iatom < _nbeads; ++iatom) {
            Force = conf->getBead(iatom)->getF();
            _b(_least_sq_offset + 3 * _nbeads * _frame_counter + iatom) = Force.x();
            _b(_least_sq_offset + 3 * _nbeads * _frame_counter + _nbeads + iatom) = Force.y();
            _b(_least_sq_offset + 3 * _nbeads * _frame_counter + 2 * _nbeads + iatom) = Force.z();
        }
    } else {
        cerr << "\nERROR in csg_fmatch::EvalConfiguration - No forces in configuration!\n" << endl;
        exit(-1);
    }
    // update the frame counter
    _frame_counter += 1; 

    if (_frame_counter % _nframes == 0) { // at this point we processed _nframes frames, which is enough for one block
        // update block counter
        _nblocks++;
        // solve FM equations and accumulate the result
        FmatchAccumulateData();
        // print status information
        cout << "\nBlock No" << _nblocks << " done!" << endl;
        // write results to output files
        WriteOutFiles();

        // we must count frames from zero again for the next block
        _frame_counter = 0;
        if (_constr_least_sq) { //Constrained Least Squares
            // Matrices should be cleaned after each block is evaluated
            _A.setZero();
            _b.setZero();
            // clear and assign smoothing conditions to _B_constr
            FmatchAssignSmoothCondsToMatrix(_B_constr);
        } else { // Simple Least Squares
            // Matrices should be cleaned after each block is evaluated            
            // clear and assign smoothing conditions to _A
            FmatchAssignSmoothCondsToMatrix(_A);
            _b.setZero();
        }
    }
    if(_has_existing_forces)
        _trjreader_force->NextFrame(_top_force);
}

void CGForceMatching::FmatchAccumulateData() 
{
    if (_constr_least_sq) { // Constrained Least Squares
        // Solving linear equations system
        votca::tools::linalg_constrained_qrsolve(_x, _A, _b, _B_constr);
    } else { // Simple Least Squares
        
        Eigen::HouseholderQR<Eigen::MatrixXd> dec(_A);
        _x=dec.solve(_b);
        Eigen::VectorXd residual=_b-_A*_x;
        // calculate FM residual - quality of FM
        // FM residual is initially calculated in (kJ/(mol*nm))^2
        double fm_resid = residual.cwiseAbs2().sum();

        // strange number is units conversion -> now (kcal/(mol*angstrom))^2
        fm_resid /= 3 * _nbeads * _frame_counter * 1750.5856;

        cout << endl;
        cout << "#### Force matching residual ####" << endl;
        cout << "     Chi_2 = " << fm_resid << endl;
        cout << "#################################" << endl;
        cout << endl;
    }

    SplineContainer::iterator is;
    for (is = _splines.begin(); is != _splines.end(); ++is) {
        int &mp = (*is)->matr_pos;
        int &ngp = (*is)->num_gridpoints;

        // _x contains results for all splines. Here we cut the results for one spline
        for (int i = 0; i < ngp; i++) {
            (*is)->block_res_f[i] = _x[ i + mp ];
            (*is)->block_res_f2[i] = _x[ i + mp + ngp];
        }
        // result cutted before is assigned to the corresponding spline
        (*is)->Spline.setSplineData((*is)->block_res_f, (*is)->block_res_f2);

        // first output point = first grid point
        double out_x = (*is)->Spline.getGridPoint(0);

        // point in the middle of the output grid for printing debug information
        int grid_point_debug = (*is)->num_outgrid / 2;

        // loop over output grid
        for (int i = 0; i < (*is)->num_outgrid; i++) {
            // update resSum (add result of a particular block)
            (*is)->resSum[i] += (*is)->Spline.Calculate(out_x);
            // print useful debug information
            if (i == grid_point_debug) cout << "This should be a number: " << (*is)->Spline.Calculate(out_x) << " " << endl;
            // update resSum2 (add result of a particular block)
            (*is)->resSum2[i] += (*is)->Spline.Calculate(out_x) * (*is)->Spline.Calculate(out_x);
            // output point for the next iteration
            out_x += (*is)->dx_out;
        }
    }
}

void CGForceMatching::FmatchAssignSmoothCondsToMatrix(Eigen::MatrixXd &Matrix)
{
// This function assigns Spline smoothing conditions to the Matrix.
// For the simple least squares the function is used for matrix _A
// For constrained least squares - for matrix _B_constr
    int line_tmp, col_tmp;
    line_tmp = 0;
    col_tmp = 0;

    Matrix.setZero();


    SplineContainer::iterator is;
    for (is = _splines.begin(); is != _splines.end(); ++is) {
        int sfnum = (*is)->num_splinefun;
        (*is)->Spline.AddBCToFitMatrix(Matrix, line_tmp, col_tmp);
        //if periodic potential, one additional constraint has to be taken into account!
        if ((*is)->periodic != 0){
            (*is)->Spline.AddBCSumZeroToFitMatrix(Matrix, line_tmp, col_tmp);
            //update counter
            line_tmp += 1;
        }
        // update counters
        line_tmp += sfnum + 1;
        col_tmp += 2 * (sfnum + 1);
    }
}

void CGForceMatching::LoadOptions(const string &file) 
{
    load_property_from_xml(_options, file);
    _bonded = _options.Select("cg.bonded");
    _nonbonded = _options.Select("cg.non-bonded");
}

void CGForceMatching::EvalBonded(Topology *conf, SplineInfo *sinfo) 
{
    std::list<Interaction *> interList;
    std::list<Interaction *>::iterator interListIter;

    interList = conf->InteractionsInGroup(sinfo->splineName);

    for (interListIter = interList.begin(); interListIter != interList.end(); ++interListIter) {

        int beads_in_int = (*interListIter)->BeadCount(); // 2 for bonds, 3 for angles, 4 for dihedrals

        CubicSpline &SP = sinfo->Spline;

        int &mpos = sinfo->matr_pos;

        double var = (*interListIter)->EvaluateVar(*conf); // value of bond, angle, or dihedral

        for (int loop = 0; loop < beads_in_int; loop++) {
            int ii = (*interListIter)->getBeadId(loop);
            vec gradient = (*interListIter)->Grad(*conf, loop);

            SP.AddToFitMatrix(_A, var,
                    _least_sq_offset + 3 * _nbeads * _frame_counter + ii, mpos, -gradient.x());
            SP.AddToFitMatrix(_A, var,
                    _least_sq_offset + 3 * _nbeads * _frame_counter + _nbeads + ii, mpos, -gradient.y());
            SP.AddToFitMatrix(_A, var,
                    _least_sq_offset + 3 * _nbeads * _frame_counter + 2 * _nbeads + ii, mpos, -gradient.z());
        }
    }
}

void CGForceMatching::EvalNonbonded(Topology *conf, SplineInfo *sinfo) 
{
    // generate the neighbour list
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

   nb->setCutoff(sinfo->_options->get("fmatch.max").as<double>()); // implement different cutoffs for different interactions!

    // generate the bead lists
    BeadList beads1, beads2;
    beads1.Generate(*conf, sinfo->type1);
    beads2.Generate(*conf, sinfo->type2);

    // is it same types or different types?
    if (sinfo->type1 == sinfo->type2)
        nb->Generate(beads1, true);
    else
        nb->Generate(beads1, beads2, true);

    NBList::iterator pair_iter;
    // iterate over all pairs
    for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {
        int iatom = (*pair_iter)->first->getId();
        int jatom = (*pair_iter)->second->getId();
        double var = (*pair_iter)->dist();
        vec gradient = (*pair_iter)->r();
        gradient.normalize();

        CubicSpline &SP = sinfo->Spline;

        int &mpos = sinfo->matr_pos;

        // add iatom
        SP.AddToFitMatrix(_A, var,
                _least_sq_offset + 3 * _nbeads * _frame_counter + iatom, mpos, gradient.x());
        SP.AddToFitMatrix(_A, var,
                _least_sq_offset + 3 * _nbeads * _frame_counter + _nbeads + iatom, mpos, gradient.y());
        SP.AddToFitMatrix(_A, var,
                _least_sq_offset + 3 * _nbeads * _frame_counter + 2 * _nbeads + iatom, mpos, gradient.z());

        // add jatom 
        SP.AddToFitMatrix(_A, var,
                _least_sq_offset + 3 * _nbeads * _frame_counter + jatom, mpos, -gradient.x());
        SP.AddToFitMatrix(_A, var,
                _least_sq_offset + 3 * _nbeads * _frame_counter + _nbeads + jatom, mpos, -gradient.y());
        SP.AddToFitMatrix(_A, var,
                _least_sq_offset + 3 * _nbeads * _frame_counter + 2 * _nbeads + jatom, mpos, -gradient.z());
    }
    delete nb;
}
