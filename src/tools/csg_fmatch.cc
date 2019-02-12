/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include "csg_fmatch.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/nblistgrid.h>
#include <votca/csg/nblistgrid_3body.h>
#include <votca/tools/cubicspline.h>
#include <votca/tools/linalg.h>
#include <votca/tools/table.h>

int main(int argc, char **argv) {
  CGForceMatching app;
  return app.Exec(argc, argv);
}

void CGForceMatching::Initialize(void) {
  CsgApplication::Initialize();
  AddProgramOptions()("options", boost::program_options::value<string>(),
                      "  options file for coarse graining")(
      "trj-force", boost::program_options::value<string>(),
      "  coarse-grained trajectory containing forces of already known "
      "interactions");
}

bool CGForceMatching::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  CheckRequired("trj", "no trajectory file specified");
  CheckRequired("options", "need to specify options file");
  LoadOptions(OptionsMap()["options"].as<string>());

  _has_existing_forces = false;
  if (OptionsMap().count("trj-force"))
    _has_existing_forces = true;
  return true;
}

void CGForceMatching::BeginEvaluate(Topology *top, Topology *top_atom) {
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

  // read _nframes from input file
  _nframes = _options.get("cg.fmatch.frames_per_block").as<int>();
  // read _constr_least_sq from input file
  _constr_least_sq = _options.get("cg.fmatch.constrainedLS").as<bool>();

  // initializing bonded interactions
  for (list<Property *>::iterator iter = _bonded.begin(); iter != _bonded.end();
       ++iter) {
    SplineInfo *i = new SplineInfo(_splines.size(), true, _col_cntr, *iter);
    // adjust initial matrix dimensions:
    _line_cntr += i->num_gridpoints;
    _col_cntr += 2 * i->num_gridpoints;
    // if periodic potential, one additional constraint has to be taken into
    // account -> 1 additional line in matrix
    if (i->periodic != 0) {
      _line_cntr += 1;
    }
    // add spline to container
    _splines.push_back(i);
  }

  // initializing non-bonded interactions
  for (list<Property *>::iterator iter = _nonbonded.begin();
       iter != _nonbonded.end(); ++iter) {
    SplineInfo *i = new SplineInfo(_splines.size(), false, _col_cntr, *iter);
    // adjust initial matrix dimensions:
    // number of constraints/restraints
    _line_cntr += i->num_gridpoints;
    // number of coefficients
    _col_cntr += 2 * i->num_gridpoints;

    // preliminary: use also spline functions for the threebody interaction. So
    // far only angular interaction implemented
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
    _B_constr = Eigen::MatrixXd::Zero(_line_cntr, _col_cntr);

    // resize matrix _A
    _A = Eigen::MatrixXd::Zero(3 * _nbeads * _nframes, _col_cntr);
    // resize vector _b
    _b = Eigen::VectorXd::Zero(3 * _nbeads * _nframes);

    // in case of constrained least squares smoothing conditions
    // are assigned to matrix _B_constr
    FmatchAssignSmoothCondsToMatrix(_B_constr);
  } else { // Simple Least Squares

    cout << "\nUsing simple Least Squares! " << endl;
    // assign _least_sq_offset
    _least_sq_offset = _line_cntr;

    // resize matrix _A
    _A = Eigen::MatrixXd::Zero(_line_cntr + 3 * _nbeads * _nframes, _col_cntr);
    // resize vector _b
    _b = Eigen::VectorXd::Zero(_line_cntr + 3 * _nbeads * _nframes);

    // in case of simple least squares smoothing conditions
    // are assigned to matrix _A
    FmatchAssignSmoothCondsToMatrix(_A);
    // clear _b (only necessary in simple least squares)
    _b.setZero();
  }
  // resize and clear _x
  _x = Eigen::VectorXd::Zero(_col_cntr);

  if (_has_existing_forces) {
    _top_force.CopyTopologyData(top);
    _trjreader_force =
        TrjReaderFactory().Create(_op_vm["trj-force"].as<string>());
    if (_trjreader_force == NULL)
      throw runtime_error(string("input format not supported: ") +
                          _op_vm["trj-force"].as<string>());
    // open the trajectory
    _trjreader_force->Open(_op_vm["trj-force"].as<string>());
    // read in first frame
    _trjreader_force->FirstFrame(_top_force);
  }
}

CGForceMatching::SplineInfo::SplineInfo(int index, bool bonded_, int matr_pos_,
                                        Property *options) {
  // initialize standard data
  splineIndex = index;
  _options = options;
  splineName = options->get("name").value();
  bonded = bonded_;
  // in general natural boundary conditions are used for splines (default is no)
  periodic = 0;
  // check if non-bonded 3-body interaction or not (default is no)
  threebody = 0;
  // initialize additional parameters for threebody interactions
  //(values of Molinero water potential)
  a = 0.37;     //(0.37 nm)
  sigma = 1;    //(dimensionless)
  gamma = 0.12; //(0.12 nm = 1.2 Ang)

  // get non-bonded information
  if (!bonded) {
    // check if option threebody exists
    if (options->exists("threebody")) {
      threebody = options->get("threebody").as<bool>();
    }
    // check if threebody interaction or not
    if (threebody) {
      type1 = options->get("type1").value();
      type2 = options->get("type2").value();
      type3 = options->get("type3").value();
      // read in additional parameters for threebody interactions
      if (options->exists("fmatch.a")) {
        a = options->get("fmatch.a").as<double>();
      }
      if (options->exists("fmatch.sigma")) {
        sigma = options->get("fmatch.sigma").as<double>();
      }
      if (options->exists("fmatch.gamma")) {
        gamma = options->get("fmatch.gamma").as<double>();
      }
    }
    // if not threebody only read in the two bead types
    if (!threebody) {
      type1 = options->get("type1").value();
      type2 = options->get("type2").value();
    }
  }
  if (bonded) {
    // check if option periodic exists
    if (options->exists("fmatch.periodic")) {
      periodic = options->get("fmatch.periodic").as<bool>();
      // set cubic spline Spline boundary conditions to periodic
      Spline.setBCInt(1);
    }
  }
  std::cout << "a: " << a << " ,sigma: " << sigma << " ,gamma: " << gamma
            << std::endl;

  // initialize the grid
  double grid_min = options->get("fmatch.min").as<double>();
  double grid_max = options->get("fmatch.max").as<double>();
  double grid_step = options->get("fmatch.step").as<double>();

  // GenerateGrid returns number of grid points. We subtract 1 to get
  // the number of spline functions
  num_gridpoints = Spline.GenerateGrid(grid_min, grid_max, grid_step);
  num_splinefun = num_gridpoints - 1;

  cout << "Number of spline functions for the interaction " << splineName << ":"
       << num_splinefun << endl;

  matr_pos = matr_pos_;

  // initialize grid for block averaging
  dx_out = options->get("fmatch.out_step").as<double>();
  // number of output grid points
  num_outgrid = 1 + (int)((grid_max - grid_min) / dx_out);
  result = Eigen::VectorXd::Zero(num_outgrid);
  error = Eigen::VectorXd::Zero(num_outgrid);
  resSum = Eigen::VectorXd::Zero(num_outgrid);
  resSum2 = Eigen::VectorXd::Zero(num_outgrid);

  // Only if threebody interaction, the derivatives are explicitly calculated
  if (threebody) {
    resultDer = Eigen::VectorXd::Zero(num_outgrid);
    errorDer = Eigen::VectorXd::Zero(num_outgrid);
    resSumDer = Eigen::VectorXd::Zero(num_outgrid);
    resSumDer2 = Eigen::VectorXd::Zero(num_outgrid);
  }

  block_res_f = Eigen::VectorXd::Zero(num_outgrid);
  block_res_f2 = Eigen::VectorXd::Zero(num_outgrid);
}
void CGForceMatching::EndEvaluate() {
  // sanity check
  if (_nblocks == 0) {
    cerr << "\nERROR in CGForceMatching::EndEvaluate - No blocks have been "
            "processed so far"
         << endl;
    cerr << "It might be that you are using trajectory, which is smaller than "
            "needed for one block"
         << endl;
    cerr << "Check your input!" << endl;
    exit(-1);
  }

  cout << "\nWe are done, thank you very much!" << endl;
  if (_has_existing_forces) {
    _trjreader_force->Close();
    delete _trjreader_force;
  }
}

void CGForceMatching::WriteOutFiles() {
  string file_extension = ".force";
  string file_extension_pot = ".pot";
  string file_name;
  string file_nameDer;
  Table force_tab;
  Table force_tabDer;

  // table with error column
  force_tab.SetHasYErr(true);
  force_tabDer.SetHasYErr(true);

  SplineContainer::iterator is;

  for (is = _splines.begin(); is != _splines.end(); ++is) {
    // construct meaningful outfile name
    file_name = (*is)->splineName;

    // resize table
    force_tab.resize((*is)->num_outgrid);

    // If not threebody, the result represents the force
    if (!((*is)->threebody)) {
      file_name = file_name + file_extension;
      // print output file names on stdout
      cout << "Updating file: " << file_name << endl;
    }

    // If threebody interaction, the result represents the potential and (-1)
    // the derivative represents the force Only then, the derivatives are
    // explicitly calculated
    if ((*is)->threebody) {
      file_name = file_name + file_extension_pot;
      file_nameDer = (*is)->splineName;
      file_nameDer = file_nameDer + file_extension;

      force_tabDer.resize((*is)->num_outgrid);
      // print output file names on stdout
      cout << "Updating files: " << file_name << " and: " << file_nameDer
           << endl;
    }

    // loop over output grid points
    for (int i = 0; i < (*is)->num_outgrid; i++) {
      // average value
      (*is)->result[i] = (*is)->resSum[i] / _nblocks;
      // standard deviation of the average
      (*is)->error[i] = sqrt(
          ((*is)->resSum2[i] / _nblocks - (*is)->result[i] * (*is)->result[i]) /
          _nblocks);

      // Only if threebody interaction, the derivatives are explicitly
      // calculated
      if ((*is)->threebody) {
        // average value
        (*is)->resultDer[i] = (*is)->resSumDer[i] / _nblocks;
        // standard deviation of the average
        (*is)->errorDer[i] = sqrt(((*is)->resSumDer2[i] / _nblocks -
                                   (*is)->resultDer[i] * (*is)->resultDer[i]) /
                                  _nblocks);
      }
    }

    // first output point = first grid point
    double out_x = (*is)->Spline.getGridPoint(0);
    // loop over output grid
    for (int i = 0; i < (*is)->num_outgrid; i++) {

      // If not threebody the result is (-1) the force
      if (!((*is)->threebody)) {
        // put point, result, flag and accuracy at point out_x into the table
        force_tab.set(i, out_x, (-1.0) * (*is)->result[i], 'i',
                      (*is)->error[i]);
      }

      // If threebody interaction, force_tab represents the potential (-1) which
      // is the Antiderivative of the force Only if threebody interaction, the
      // derivatives are explicitly calculated
      if ((*is)->threebody) {
        // put point, result, flag and accuracy at point out_x into the table
        force_tab.set(i, out_x, (+1.0) * (*is)->result[i], 'i',
                      (*is)->error[i]);
        force_tabDer.set(i, out_x, (-1.0) * (*is)->resultDer[i], 'i',
                         (*is)->errorDer[i]);
      }

      // update out_x for the next iteration
      out_x += (*is)->dx_out;
    }
    // save table in the file
    force_tab.Save(file_name);

    // clear the table for the next spline
    force_tab.clear();

    // Only if threebody interaction, the derivatives are explicitly calculated
    if ((*is)->threebody) {
      force_tabDer.Save(file_nameDer);
      // clear the table for the next spline
      force_tabDer.clear();
    }
  }
}

void CGForceMatching::EvalConfiguration(Topology *conf, Topology *conf_atom) {
  SplineContainer::iterator spiter;
  if (conf->BeadCount() == 0)
    throw std::runtime_error(
        "CG Topology has 0 beads, check your mapping file!");
  if (_has_existing_forces) {
    if (conf->BeadCount() != _top_force.BeadCount())
      throw std::runtime_error(
          "number of beads in topology and force topology does not match");
    for (int i = 0; i < conf->BeadCount(); ++i) {
      conf->getBead(i)->F() -= _top_force.getBead(i)->getF();
      vec d = conf->getBead(i)->getPos() - _top_force.getBead(i)->getPos();
      if (abs(d) > _dist) { // default is 1e-5, otherwise it can be a too
                            // strict criterion
        throw std::runtime_error(
            "One or more bead positions in mapped and reference force "
            "trajectory differ by more than 1e-5");
      }
    }
  }

  for (spiter = _splines.begin(); spiter != _splines.end(); ++spiter) {
    SplineInfo *sinfo = *spiter;
    if (sinfo->bonded) // bonded interaction
      EvalBonded(conf, sinfo);
    else // non-bonded interaction
        // check if threebody interaction or not
        if (sinfo->threebody) {
      EvalNonbonded_Threebody(conf, sinfo);
    } else {
      EvalNonbonded(conf, sinfo);
    }
  }

  // loop for the forces vector:
  // hack, chage the Has functions..
  if (conf->getBead(0)->HasF()) {
    vec Force(0., 0., 0.);
    for (int iatom = 0; iatom < _nbeads; ++iatom) {
      Force = conf->getBead(iatom)->getF();
      _b(_least_sq_offset + 3 * _nbeads * _frame_counter + iatom) = Force.x();
      _b(_least_sq_offset + 3 * _nbeads * _frame_counter + _nbeads + iatom) =
          Force.y();
      _b(_least_sq_offset + 3 * _nbeads * _frame_counter + 2 * _nbeads +
         iatom) = Force.z();
    }
  } else {
    cerr << "\nERROR in csg_fmatch::EvalConfiguration - No forces in "
            "configuration!\n"
         << endl;
    exit(-1);
  }
  // update the frame counter
  _frame_counter += 1;

  if (_frame_counter % _nframes == 0) { // at this point we processed _nframes
                                        // frames, which is enough for one
                                        // block
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
    if (_constr_least_sq) { // Constrained Least Squares
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
  if (_has_existing_forces)
    _trjreader_force->NextFrame(_top_force);
}

void CGForceMatching::FmatchAccumulateData() {
  if (_constr_least_sq) { // Constrained Least Squares
    // Solving linear equations system
    votca::tools::linalg_constrained_qrsolve(_x, _A, _b, _B_constr);
  } else { // Simple Least Squares

    Eigen::HouseholderQR<Eigen::MatrixXd> dec(_A);
    _x = dec.solve(_b);
    Eigen::VectorXd residual = _b - _A * _x;
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

    // _x contains results for all splines. Here we cut the results for one
    // spline
    for (int i = 0; i < ngp; i++) {
      (*is)->block_res_f[i] = _x[i + mp];
      (*is)->block_res_f2[i] = _x[i + mp + ngp];
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
      // update resSum2 (add result of a particular block)
      (*is)->resSum2[i] +=
          (*is)->Spline.Calculate(out_x) * (*is)->Spline.Calculate(out_x);

      // Only if threebody interaction, the derivatives are explicitly
      // calculated
      if ((*is)->threebody) {
        (*is)->resSumDer[i] += (*is)->Spline.CalculateDerivative(out_x);
        // update resSumDer2 (add result of a particular block)
        (*is)->resSumDer2[i] += (*is)->Spline.CalculateDerivative(out_x) *
                                (*is)->Spline.CalculateDerivative(out_x);
      }

      // print useful debug information
      if (i == grid_point_debug)
        cout << "This should be a number: " << (*is)->Spline.Calculate(out_x)
             << " " << endl;

      // output point for the next iteration
      out_x += (*is)->dx_out;
    }
  }
}

void CGForceMatching::FmatchAssignSmoothCondsToMatrix(Eigen::MatrixXd &Matrix) {
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
    // if periodic potential, one additional constraint has to be taken into
    // account!
    if ((*is)->periodic != 0) {
      (*is)->Spline.AddBCSumZeroToFitMatrix(Matrix, line_tmp, col_tmp);
      // update counter
      line_tmp += 1;
    }

    // update counters
    line_tmp += sfnum + 1;
    col_tmp += 2 * (sfnum + 1);
  }
}

void CGForceMatching::LoadOptions(const string &file) {
  load_property_from_xml(_options, file);
  _bonded = _options.Select("cg.bonded");
  _nonbonded = _options.Select("cg.non-bonded");
}

void CGForceMatching::EvalBonded(Topology *conf, SplineInfo *sinfo) {
  std::list<Interaction *> interList;
  std::list<Interaction *>::iterator interListIter;

  interList = conf->InteractionsInGroup(sinfo->splineName);

  for (interListIter = interList.begin(); interListIter != interList.end();
       ++interListIter) {

    int beads_in_int =
        (*interListIter)->BeadCount(); // 2 for bonds, 3 for angles, 4 for
                                       // dihedrals

    CubicSpline &SP = sinfo->Spline;

    int &mpos = sinfo->matr_pos;

    double var = (*interListIter)->EvaluateVar(*conf); // value of bond, angle,
                                                       // or dihedral

    for (int loop = 0; loop < beads_in_int; loop++) {
      int ii = (*interListIter)->getBeadId(loop);
      vec gradient = (*interListIter)->Grad(*conf, loop);

      SP.AddToFitMatrix(_A, var,
                        _least_sq_offset + 3 * _nbeads * _frame_counter + ii,
                        mpos, -gradient.x());
      SP.AddToFitMatrix(_A, var,
                        _least_sq_offset + 3 * _nbeads * _frame_counter +
                            _nbeads + ii,
                        mpos, -gradient.y());
      SP.AddToFitMatrix(_A, var,
                        _least_sq_offset + 3 * _nbeads * _frame_counter +
                            2 * _nbeads + ii,
                        mpos, -gradient.z());
    }
  }
}

void CGForceMatching::EvalNonbonded(Topology *conf, SplineInfo *sinfo) {

  // generate the neighbour list
  NBList *nb;

  bool gridsearch = false;

  if (_options.exists("cg.nbsearch")) {
    if (_options.get("cg.nbsearch").as<string>() == "grid")
      gridsearch = true;
    else if (_options.get("cg.nbsearch").as<string>() == "simple")
      gridsearch = false;
    else
      throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
  }
  if (gridsearch)
    nb = new NBListGrid();
  else
    nb = new NBList();

  nb->setCutoff(
      sinfo->_options->get("fmatch.max").as<double>()); // implement different
                                                        // cutoffs for
                                                        // different
                                                        // interactions!

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
                      _least_sq_offset + 3 * _nbeads * _frame_counter + iatom,
                      mpos, gradient.x());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   _nbeads + iatom,
                      mpos, gradient.y());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   2 * _nbeads + iatom,
                      mpos, gradient.z());

    // add jatom
    SP.AddToFitMatrix(_A, var,
                      _least_sq_offset + 3 * _nbeads * _frame_counter + jatom,
                      mpos, -gradient.x());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   _nbeads + jatom,
                      mpos, -gradient.y());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   2 * _nbeads + jatom,
                      mpos, -gradient.z());
  }
  delete nb;
}

void CGForceMatching::EvalNonbonded_Threebody(Topology *conf,
                                              SplineInfo *sinfo) {
  // so far option gridsearch ignored. Only simple search

  // generate the neighbour list
  NBList_3Body *nb;

  bool gridsearch = false;

  if (_options.exists("cg.nbsearch")) {
    if (_options.get("cg.nbsearch").as<string>() == "grid")
      gridsearch = true;
    else if (_options.get("cg.nbsearch").as<string>() == "simple")
      gridsearch = false;
    else
      throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
  }
  if (gridsearch)
    nb = new NBListGrid_3Body();
  else
    nb = new NBList_3Body();

  nb->setCutoff(sinfo->a); // implement different cutoffs for different
                           // interactions!
  // Here, a is the distance between two beads of a triple, where the 3-body
  // interaction is zero

  // generate the bead lists
  BeadList beads1, beads2, beads3;
  beads1.Generate(*conf, sinfo->type1);
  beads2.Generate(*conf, sinfo->type2);
  beads3.Generate(*conf, sinfo->type3);

  // check if type1 and type2 are the same
  if (sinfo->type1 == sinfo->type2) {
    // if all three types are the same
    if (sinfo->type2 == sinfo->type3) {
      nb->Generate(beads1, true);
    }
    // if type2 and type3 are different, use the Generate function for 2 bead
    // types
    if (sinfo->type2 != sinfo->type3) {
      nb->Generate(beads1, beads3, true);
    }
  }
  // if type1 != type2
  if (sinfo->type1 != sinfo->type2) {
    // if the last two types are the same, use Generate function with them as
    // the first two bead types Neighborlist_3body is constructed in a way that
    // the two equal bead types have two be the first 2 types
    if (sinfo->type2 == sinfo->type3) {
      nb->Generate(beads1, beads2, true);
    }
    if (sinfo->type2 != sinfo->type3) {
      // type1 = type3 !=type2
      if (sinfo->type1 == sinfo->type3) {
        nb->Generate(beads2, beads1, true);
      }
      // type1 != type2 != type3
      if (sinfo->type1 != sinfo->type3) {
        nb->Generate(beads1, beads2, beads3, true);
      }
    }
  }

  NBList_3Body::iterator triple_iter;
  // iterate over all triples
  for (triple_iter = nb->begin(); triple_iter != nb->end(); ++triple_iter) {
    int iatom = (*triple_iter)->bead1()->getId();
    int jatom = (*triple_iter)->bead2()->getId();
    int katom = (*triple_iter)->bead3()->getId();
    double distij = (*triple_iter)->dist12();
    double distik = (*triple_iter)->dist13();
    vec rij = (*triple_iter)->r12();
    vec rik = (*triple_iter)->r13();

    double gamma_sigma = (sinfo->gamma) * (sinfo->sigma);
    double denomij = (distij - (sinfo->a) * (sinfo->sigma));
    double denomik = (distik - (sinfo->a) * (sinfo->sigma));
    double expij = exp(gamma_sigma / denomij);
    double expik = exp(gamma_sigma / denomik);

    vec gradient1, gradient2;

    CubicSpline &SP = sinfo->Spline;

    int &mpos = sinfo->matr_pos;

    double var = acos(rij * rik / sqrt((rij * rij) * (rik * rik)));

    double acos_prime =
        1.0 /
        (sqrt(1 -
              (rij * rik) * (rij * rik) / (distij * distik * distij * distik)));

    // evaluate gradient1 and gradient2 for iatom:
    gradient1 = acos_prime *
                ((rij + rik) / (distij * distik) -
                 (rij * rik) * ((rik * rik) * rij + (rij * rij) * rik) /
                     (distij * distij * distij * distik * distik * distik)) *
                expij * expik;
    // gradient2
    gradient2 = ((rij / distij) * (gamma_sigma / (denomij * denomij)) +
                 (rik / distik) * (gamma_sigma / (denomik * denomik))) *
                expij * expik;

    // add iatom
    SP.AddToFitMatrix(_A, var,
                      _least_sq_offset + 3 * _nbeads * _frame_counter + iatom,
                      mpos, -gradient1.x(), -gradient2.x());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   _nbeads + iatom,
                      mpos, -gradient1.y(), -gradient2.y());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   2 * _nbeads + iatom,
                      mpos, -gradient1.z(), -gradient2.z());

    // evaluate gradient1 and gradient2 for jatom:
    gradient1 = acos_prime *
                (-rik / (distij * distik) +
                 (rij * rik) * rij / (distik * distij * distij * distij)) *
                expij * expik;
    // gradient2
    gradient2 = ((rij / distij) * (-1.0 * gamma_sigma / (denomij * denomij))) *
                expij * expik;

    // add jatom
    SP.AddToFitMatrix(_A, var,
                      _least_sq_offset + 3 * _nbeads * _frame_counter + jatom,
                      mpos, -gradient1.x(), -gradient2.x());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   _nbeads + jatom,
                      mpos, -gradient1.y(), -gradient2.y());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   2 * _nbeads + jatom,
                      mpos, -gradient1.z(), -gradient2.z());

    // evaluate gradient1 and gradient2 for katom:
    gradient1 = acos_prime *
                (-rij / (distij * distik) +
                 (rij * rik) * rik / (distij * distik * distik * distik)) *
                expij * expik;
    // gradient2
    gradient2 = ((rik / distik) * (-1.0 * gamma_sigma / (denomik * denomik))) *
                expij * expik;

    // add jatom
    SP.AddToFitMatrix(_A, var,
                      _least_sq_offset + 3 * _nbeads * _frame_counter + katom,
                      mpos, -gradient1.x(), -gradient2.x());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   _nbeads + katom,
                      mpos, -gradient1.y(), -gradient2.y());
    SP.AddToFitMatrix(_A, var, _least_sq_offset + 3 * _nbeads * _frame_counter +
                                   2 * _nbeads + katom,
                      mpos, -gradient1.z(), -gradient2.z());
  }
  delete nb;
}
