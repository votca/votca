/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

// VOTCA includes
#include <votca/tools/cubicspline.h>
#include <votca/tools/linalg.h>
#include <votca/tools/table.h>

// Local VOTCA includes
#include "votca/csg/beadlist.h"
#include "votca/csg/nblistgrid.h"
#include "votca/csg/nblistgrid_3body.h"

// Local private includes
#include "csg_fmatch.h"

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

  has_existing_forces_ = false;
  if (OptionsMap().count("trj-force")) {
    has_existing_forces_ = true;
  }
  return true;
}

void CGForceMatching::BeginEvaluate(Topology *top, Topology *) {
  // set counters to zero value:
  nblocks_ = 0;
  line_cntr_ = col_cntr_ = 0;

  // Number of CG beads in topology
  nbeads_ = top->BeadCount();
  // Set frame counter to zero
  frame_counter_ = 0;

  // accuracy for evaluating the difference in bead positions (default 1e-5)
  dist_ = 1e-5;
  if (options_.exists("cg.fmatch.dist")) {
    dist_ = options_.get("cg.fmatch.dist").as<double>();
  }

  // read  nframes_ from input file
  nframes_ = options_.get("cg.fmatch.frames_per_block").as<votca::Index>();
  // read  constr_least_sq_ from input file
  constr_least_sq_ = options_.get("cg.fmatch.constrainedLS").as<bool>();

  // initializing bonded interactions
  for (votca::tools::Property *prop : bonded_) {
    // add spline to container
    splines_.emplace_back(splines_.size(), true, col_cntr_, prop);
    // adjust initial Eigen::Matrix3d dimensions:
    line_cntr_ += splines_.back().num_gridpoints;
    col_cntr_ += 2 * splines_.back().num_gridpoints;
    // if periodic potential, one additional constraint has to be taken into
    // account -> 1 additional line in matrix
    if (splines_.back().periodic != 0) {
      line_cntr_ += 1;
    }
  }

  // initializing non-bonded interactions
  for (votca::tools::Property *prop : nonbonded_) {
    // add spline to container
    splines_.emplace_back(splines_.size(), false, col_cntr_, prop);
    // adjust initial Eigen::Matrix3d dimensions:
    // number of constraints/restraints
    line_cntr_ += splines_.back().num_gridpoints;
    // number of coefficients
    col_cntr_ += 2 * splines_.back().num_gridpoints;

    // preliminary: use also spline functions for the threebody interaction. So
    // far only angular interaction implemented
  }

  cout << "\nYou are using VOTCA!\n";
  cout << "\nhey, somebody wants to forcematch!\n";

  // now initialize  A_,  b_,  x_ and probably  B_constr_
  // depending on least-squares algorithm used
  if (constr_least_sq_) {  // Constrained Least Squares

    cout << "\nUsing constrained Least Squares!\n " << endl;

    // assign  least_sq_offset_
    least_sq_offset_ = 0;

    // resize and clear  B_constr_
    B_constr_ = Eigen::MatrixXd::Zero(line_cntr_, col_cntr_);

    // resize Eigen::Matrix3d  A_
    A_ = Eigen::MatrixXd::Zero(3 * nbeads_ * nframes_, col_cntr_);
    // resize vector  b_
    b_ = Eigen::VectorXd::Zero(3 * nbeads_ * nframes_);

    // in case of constrained least squares smoothing conditions
    // are assigned to Eigen::Matrix3d  B_constr_
    FmatchAssignSmoothCondsToMatrix(B_constr_);
  } else {  // Simple Least Squares

    cout << "\nUsing simple Least Squares! " << endl;
    // assign  least_sq_offset_
    least_sq_offset_ = line_cntr_;

    // resize Eigen::Matrix3d  A_
    A_ = Eigen::MatrixXd::Zero(line_cntr_ + 3 * nbeads_ * nframes_, col_cntr_);
    // resize vector  b_
    b_ = Eigen::VectorXd::Zero(line_cntr_ + 3 * nbeads_ * nframes_);

    // in case of simple least squares smoothing conditions
    // are assigned to Eigen::Matrix3d  A_
    FmatchAssignSmoothCondsToMatrix(A_);
    // clear  b_ (only necessary in simple least squares)
    b_.setZero();
  }
  // resize and clear  x_
  x_ = Eigen::VectorXd::Zero(col_cntr_);

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

CGForceMatching::SplineInfo::SplineInfo(votca::Index index,
                                        bool bonded_interaction,
                                        votca::Index matr_pos_col,
                                        votca::tools::Property *options) {
  // initialize standard data
  splineIndex = index;
  options_ = options;
  splineName = options->get("name").value();
  bonded = bonded_interaction;
  // in general natural boundary conditions are used for splines (default is no)
  periodic = 0;
  // check if non-bonded 3-body interaction or not (default is no)
  threebody = 0;
  // initialize additional parameters for threebody interactions
  //(values of Molinero water potential)
  a = 0.37;      //(0.37 nm)
  sigma = 1;     //(dimensionless)
  gamma = 0.12;  //(0.12 nm = 1.2 Ang)

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

  matr_pos = matr_pos_col;

  // initialize grid for block averaging
  dx_out = options->get("fmatch.out_step").as<double>();
  // number of output grid points
  num_outgrid = 1 + (votca::Index)((grid_max - grid_min) / dx_out);
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
  if (nblocks_ == 0) {
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
  if (has_existing_forces_) {
    trjreader_force_->Close();
  }
}

void CGForceMatching::WriteOutFiles() {
  string file_extension = ".force";
  string file_extension_pot = ".pot";
  string file_name;
  string file_nameDer;
  votca::tools::Table force_tab;
  votca::tools::Table force_tabDer;

  // table with error column
  force_tab.SetHasYErr(true);
  force_tabDer.SetHasYErr(true);

  for (SplineInfo &spline : splines_) {
    // construct meaningful outfile name
    file_name = spline.splineName;

    // resize table
    force_tab.resize(spline.num_outgrid);

    // If not threebody, the result represents the force
    if (!(spline.threebody)) {
      file_name = file_name + file_extension;
      // print output file names on stdout
      cout << "Updating file: " << file_name << endl;
    }

    // If threebody interaction, the result represents the potential and (-1)
    // the derivative represents the force Only then, the derivatives are
    // explicitly calculated
    if (spline.threebody) {
      file_name = file_name + file_extension_pot;
      file_nameDer = spline.splineName;
      file_nameDer = file_nameDer + file_extension;

      force_tabDer.resize(spline.num_outgrid);
      // print output file names on stdout
      cout << "Updating files: " << file_name << " and: " << file_nameDer
           << endl;
    }

    spline.result = (spline.resSum).array() / nblocks_;
    spline.error = (((spline.resSum2).array() / nblocks_ -
                     (spline.result).array().abs2()) /
                    nblocks_)
                       .abs()
                       .sqrt();

    if (spline.threebody) {
      spline.resultDer = (spline.resSumDer).array() / nblocks_;
      spline.errorDer = (((spline.resSumDer2).array() / nblocks_ -
                          (spline.resultDer).array().abs2()) /
                         nblocks_)
                            .abs()
                            .sqrt();
    }

    // first output point = first grid point
    double out_x = spline.Spline.getGridPoint(0);
    // loop over output grid
    for (votca::Index i = 0; i < spline.num_outgrid; i++) {

      // If not threebody the result is (-1) the force
      if (!(spline.threebody)) {
        // put point, result, flag and accuracy at point out_x into the table
        force_tab.set(i, out_x, (-1.0) * spline.result[i], 'i',
                      spline.error[i]);
      }

      // If threebody interaction, force_tab represents the potential (-1) which
      // is the Antiderivative of the force Only if threebody interaction, the
      // derivatives are explicitly calculated
      if (spline.threebody) {
        // put point, result, flag and accuracy at point out_x into the table
        force_tab.set(i, out_x, (+1.0) * spline.result[i], 'i',
                      spline.error[i]);
        force_tabDer.set(i, out_x, (-1.0) * spline.resultDer[i], 'i',
                         spline.errorDer[i]);
      }

      // update out_x for the next iteration
      out_x += spline.dx_out;
    }
    // save table in the file
    force_tab.Save(file_name);

    // clear the table for the next spline
    force_tab.clear();

    // Only if threebody interaction, the derivatives are explicitly calculated
    if (spline.threebody) {
      force_tabDer.Save(file_nameDer);
      // clear the table for the next spline
      force_tabDer.clear();
    }
  }
}

void CGForceMatching::EvalConfiguration(Topology *conf, Topology *) {
  if (conf->BeadCount() == 0) {
    throw std::runtime_error(
        "CG Topology has 0 beads, check your mapping file!");
  }
  if (has_existing_forces_) {
    if (conf->BeadCount() != top_force_.BeadCount()) {
      throw std::runtime_error(
          "number of beads in topology and force topology does not match");
    }
    for (votca::Index i = 0; i < conf->BeadCount(); ++i) {
      conf->getBead(i)->F() -= top_force_.getBead(i)->getF();
      Eigen::Vector3d d =
          conf->getBead(i)->getPos() - top_force_.getBead(i)->getPos();
      if (d.norm() > dist_) {  // default is 1e-5, otherwise it can be a too
                               // strict criterion
        throw std::runtime_error(
            "One or more bead positions in mapped and reference force "
            "trajectory differ by more than 1e-5");
      }
    }
  }

  for (SplineInfo &sinfo : splines_) {
    if (sinfo.bonded) {
      EvalBonded(conf, &sinfo);
    } else {
      if (sinfo.threebody) {
        EvalNonbonded_Threebody(conf, &sinfo);
      } else {
        EvalNonbonded(conf, &sinfo);
      }
    }
  }

  // loop for the forces vector:
  // hack, change the Has functions..
  if (conf->getBead(0)->HasF()) {
    for (votca::Index iatom = 0; iatom < nbeads_; ++iatom) {
      const Eigen::Vector3d &Force = conf->getBead(iatom)->getF();
      b_(least_sq_offset_ + 3 * nbeads_ * frame_counter_ + iatom) = Force.x();
      b_(least_sq_offset_ + 3 * nbeads_ * frame_counter_ + nbeads_ + iatom) =
          Force.y();
      b_(least_sq_offset_ + 3 * nbeads_ * frame_counter_ + 2 * nbeads_ +
         iatom) = Force.z();
    }
  } else {
    throw std::runtime_error(
        "\nERROR in csg_fmatch::EvalConfiguration - No forces in "
        "configuration!");
  }
  // update the frame counter
  frame_counter_ += 1;

  if (frame_counter_ % nframes_ == 0) {  // at this point we processed  nframes_
                                         // frames, which is enough for one
                                         // block
                                         // update block counter
    nblocks_++;
    // solve FM equations and accumulate the result
    FmatchAccumulateData();
    // print status information
    cout << "\nBlock No" << nblocks_ << " done!" << endl;
    // write results to output files
    WriteOutFiles();

    // we must count frames from zero again for the next block
    frame_counter_ = 0;
    if (constr_least_sq_) {  // Constrained Least Squares
      // Matrices should be cleaned after each block is evaluated
      A_.setZero();
      b_.setZero();
      // clear and assign smoothing conditions to  B_constr_
      FmatchAssignSmoothCondsToMatrix(B_constr_);
    } else {  // Simple Least Squares
      // Matrices should be cleaned after each block is evaluated
      // clear and assign smoothing conditions to  A_
      FmatchAssignSmoothCondsToMatrix(A_);
      b_.setZero();
    }
  }
  if (has_existing_forces_) {
    trjreader_force_->NextFrame(top_force_);
  }
}

void CGForceMatching::FmatchAccumulateData() {
  if (constr_least_sq_) {  // Constrained Least Squares
                           // Solving linear equations system
    x_ = votca::tools::linalg_constrained_qrsolve(A_, b_, B_constr_);
  } else {  // Simple Least Squares

    Eigen::HouseholderQR<Eigen::MatrixXd> dec(A_);
    x_ = dec.solve(b_);
    Eigen::VectorXd residual = b_ - A_ * x_;
    // calculate FM residual - quality of FM
    // FM residual is calculated in (kJ/(mol*nm))^2
    double fm_resid = residual.cwiseAbs2().sum();

    fm_resid /= (double)(3 * nbeads_ * frame_counter_);

    cout << endl;
    cout << "#### Force matching residual ####" << endl;
    cout << "     Chi_2[(kJ/(mol*nm))^2] = " << fm_resid << endl;
    cout << "#################################" << endl;
    cout << endl;
  }

  for (SplineInfo &sinfo : splines_) {
    votca::Index mp = sinfo.matr_pos;
    votca::Index ngp = sinfo.num_gridpoints;

    //  x_ contains results for all splines. Here we cut the results for one
    // spline

    sinfo.block_res_f = x_.segment(mp, ngp);
    sinfo.block_res_f2 = x_.segment(mp + ngp, ngp);

    // result cut before is assigned to the corresponding spline
    sinfo.Spline.setSplineData(sinfo.block_res_f, sinfo.block_res_f2);

    // first output point = first grid point
    double out_x = sinfo.Spline.getGridPoint(0);

    // point in the middle of the output grid for printing debug information
    votca::Index grid_point_debug = sinfo.num_outgrid / 2;

    // loop over output grid
    for (votca::Index i = 0; i < sinfo.num_outgrid; i++) {
      // update resSum (add result of a particular block)
      sinfo.resSum[i] += sinfo.Spline.Calculate(out_x);
      // update resSum2 (add result of a particular block)
      sinfo.resSum2[i] +=
          sinfo.Spline.Calculate(out_x) * sinfo.Spline.Calculate(out_x);

      // Only if threebody interaction, the derivatives are explicitly
      // calculated
      if (sinfo.threebody) {
        sinfo.resSumDer[i] += sinfo.Spline.CalculateDerivative(out_x);
        // update resSumDer2 (add result of a particular block)
        sinfo.resSumDer2[i] += sinfo.Spline.CalculateDerivative(out_x) *
                               sinfo.Spline.CalculateDerivative(out_x);
      }

      // print useful debug information
      if (i == grid_point_debug) {
        cout << "This should be a number: " << sinfo.Spline.Calculate(out_x)
             << " " << endl;
      }

      // output point for the next iteration
      out_x += sinfo.dx_out;
    }
  }
}

void CGForceMatching::FmatchAssignSmoothCondsToMatrix(Eigen::MatrixXd &Matrix) {
  // This function assigns Spline smoothing conditions to the Matrix.
  // For the simple least squares the function is used for Eigen::Matrix3d  A_
  // For constrained least squares - for Eigen::Matrix3d  B_constr_

  Matrix.setZero();
  votca::Index line_tmp = 0;
  votca::Index col_tmp = 0;

  for (SplineInfo &sinfo : splines_) {

    sinfo.Spline.AddBCToFitMatrix(Matrix, line_tmp, col_tmp);
    // if periodic potential, one additional constraint has to be taken into
    // account!
    if (sinfo.periodic != 0) {
      sinfo.Spline.AddBCSumZeroToFitMatrix(Matrix, line_tmp, col_tmp);
      // update counter
      line_tmp += 1;
    }
    // update counters
    votca::Index sfnum = sinfo.num_splinefun;
    line_tmp += sfnum + 1;
    col_tmp += 2 * (sfnum + 1);
  }
}

void CGForceMatching::LoadOptions(const string &file) {
  options_.LoadFromXML(file);
  bonded_ = options_.Select("cg.bonded");
  nonbonded_ = options_.Select("cg.non-bonded");
}

void CGForceMatching::EvalBonded(Topology *conf, SplineInfo *sinfo) {

  std::vector<Interaction *> interList =
      conf->InteractionsInGroup(sinfo->splineName);

  for (Interaction *inter : interList) {

    votca::Index beads_in_int = inter->BeadCount();  // 2 for bonds, 3 for
                                                     // angles, 4 for dihedrals

    votca::tools::CubicSpline &SP = sinfo->Spline;

    votca::Index mpos = sinfo->matr_pos;

    double var = inter->EvaluateVar(*conf);  // value of bond, angle,
                                             // or dihedral

    for (votca::Index loop = 0; loop < beads_in_int; loop++) {
      votca::Index ii = inter->getBeadId(loop);
      Eigen::Vector3d gradient = inter->Grad(*conf, loop);

      SP.AddToFitMatrix(A_, var,
                        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + ii,
                        mpos, -gradient.x());
      SP.AddToFitMatrix(
          A_, var,
          least_sq_offset_ + 3 * nbeads_ * frame_counter_ + nbeads_ + ii, mpos,
          -gradient.y());
      SP.AddToFitMatrix(
          A_, var,
          least_sq_offset_ + 3 * nbeads_ * frame_counter_ + 2 * nbeads_ + ii,
          mpos, -gradient.z());
    }
  }
}

void CGForceMatching::EvalNonbonded(Topology *conf, SplineInfo *sinfo) {

  // generate the neighbour list
  std::unique_ptr<NBList> nb;

  bool gridsearch = false;

  if (options_.exists("cg.nbsearch")) {
    if (options_.get("cg.nbsearch").as<string>() == "grid") {
      gridsearch = true;
    } else if (options_.get("cg.nbsearch").as<string>() == "simple") {
      gridsearch = false;
    } else {
      throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
  }
  if (gridsearch) {
    nb = std::make_unique<NBListGrid>();
  } else {
    nb = std::make_unique<NBList>();
  }

  nb->setCutoff(
      sinfo->options_->get("fmatch.max").as<double>());  // implement different
                                                         // cutoffs for
                                                         // different
                                                         // interactions!

  // generate the bead lists
  BeadList beads1, beads2;
  beads1.Generate(*conf, sinfo->type1);
  beads2.Generate(*conf, sinfo->type2);

  // is it same types or different types?
  if (sinfo->type1 == sinfo->type2) {
    nb->Generate(beads1, true);
  } else {
    nb->Generate(beads1, beads2, true);
  }

  for (BeadPair *pair : *nb) {
    votca::Index iatom = pair->first()->getId();
    votca::Index jatom = pair->second()->getId();
    double var = pair->dist();
    Eigen::Vector3d gradient = pair->r();
    gradient.normalize();

    votca::tools::CubicSpline &SP = sinfo->Spline;

    votca::Index mpos = sinfo->matr_pos;

    // add iatom
    SP.AddToFitMatrix(A_, var,
                      least_sq_offset_ + 3 * nbeads_ * frame_counter_ + iatom,
                      mpos, gradient.x());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + nbeads_ + iatom, mpos,
        gradient.y());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + 2 * nbeads_ + iatom,
        mpos, gradient.z());

    // add jatom
    SP.AddToFitMatrix(A_, var,
                      least_sq_offset_ + 3 * nbeads_ * frame_counter_ + jatom,
                      mpos, -gradient.x());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + nbeads_ + jatom, mpos,
        -gradient.y());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + 2 * nbeads_ + jatom,
        mpos, -gradient.z());
  }
}

void CGForceMatching::EvalNonbonded_Threebody(Topology *conf,
                                              SplineInfo *sinfo) {
  // so far option gridsearch ignored. Only simple search

  // generate the neighbour list
  std::unique_ptr<NBList_3Body> nb;

  bool gridsearch = false;

  if (options_.exists("cg.nbsearch")) {
    if (options_.get("cg.nbsearch").as<string>() == "grid") {
      gridsearch = true;
    } else if (options_.get("cg.nbsearch").as<string>() == "simple") {
      gridsearch = false;
    } else {
      throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
  }
  if (gridsearch) {
    nb = std::make_unique<NBListGrid_3Body>();
  } else {
    nb = std::make_unique<NBList_3Body>();
  }

  nb->setCutoff(sinfo->a);  // implement different cutoffs for different
                            // interactions!
  // Here, a is the distance between two beads of a triple, where the 3-body
  // interaction is zero

  // generate the bead lists
  BeadList beads1, beads2, beads3;
  beads1.Generate(*conf, sinfo->type1);
  beads2.Generate(*conf, sinfo->type2);
  beads3.Generate(*conf, sinfo->type3);

  // Generate the 3body neighbour lists
  // check if type2 and type3 are the same
  if (sinfo->type2 == sinfo->type3) {
    // if then type2 and type1 are the same, all three types are the same
    // use the Generate function for this case
    if (sinfo->type1 == sinfo->type2) {
      nb->Generate(beads1, true);
    }
    // else use the Generate function for type2 being equal to type3 (and type1 being different)
    if (sinfo->type1 != sinfo->type2) {
      nb->Generate(beads1, beads2, true);
    }
  }
  // If type2 and type3 are not the same, use the Generate function for three different bead types
  // (Even if type1 and type2 or type1 and type3 are the same, the Generate function for two different beadtypes
  // is only applicable for the case that type2 is equal to type3
  if (sinfo->type2 != sinfo->type3) {
    nb->Generate(beads1, beads2, beads3, true);
  }

  for (BeadTriple *triple : *nb) {
    votca::Index iatom = triple->bead1()->getId();
    votca::Index jatom = triple->bead2()->getId();
    votca::Index katom = triple->bead3()->getId();
    double distij = triple->dist12();
    double distik = triple->dist13();
    Eigen::Vector3d rij = triple->r12();
    Eigen::Vector3d rik = triple->r13();

    double gamma_sigma = (sinfo->gamma) * (sinfo->sigma);
    double denomij = (distij - (sinfo->a) * (sinfo->sigma));
    double denomik = (distik - (sinfo->a) * (sinfo->sigma));
    double expij = std::exp(gamma_sigma / denomij);
    double expik = std::exp(gamma_sigma / denomik);

    votca::tools::CubicSpline &SP = sinfo->Spline;

    votca::Index mpos = sinfo->matr_pos;

    double var =
        std::acos(rij.dot(rik) / sqrt(rij.squaredNorm() * rik.squaredNorm()));

    double acos_prime =
        1.0 / (sqrt(1 - std::pow(rij.dot(rik), 2) /
                            (distij * distik * distij * distik)));

    Eigen::Vector3d gradient1 =
        acos_prime *
        ((rij + rik) / (distij * distik) -
         rij.dot(rik) * (rik.squaredNorm() * rij + rij.squaredNorm() * rik) /
             (distij * distij * distij * distik * distik * distik)) *
        expij * expik;
    Eigen::Vector3d gradient2 =
        ((rij / distij) * (gamma_sigma / (denomij * denomij)) +
         (rik / distik) * (gamma_sigma / (denomik * denomik))) *
        expij * expik;

    // add iatom
    SP.AddToFitMatrix(A_, var,
                      least_sq_offset_ + 3 * nbeads_ * frame_counter_ + iatom,
                      mpos, -gradient1.x(), -gradient2.x());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + nbeads_ + iatom, mpos,
        -gradient1.y(), -gradient2.y());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + 2 * nbeads_ + iatom,
        mpos, -gradient1.z(), -gradient2.z());

    // evaluate gradient1 and gradient2 for jatom:
    gradient1 = acos_prime *
                (-rik / (distij * distik) +
                 rij.dot(rik) * rij / (distik * distij * distij * distij)) *
                expij * expik;
    // gradient2
    gradient2 = ((rij / distij) * (-1.0 * gamma_sigma / (denomij * denomij))) *
                expij * expik;

    // add jatom
    SP.AddToFitMatrix(A_, var,
                      least_sq_offset_ + 3 * nbeads_ * frame_counter_ + jatom,
                      mpos, -gradient1.x(), -gradient2.x());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + nbeads_ + jatom, mpos,
        -gradient1.y(), -gradient2.y());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + 2 * nbeads_ + jatom,
        mpos, -gradient1.z(), -gradient2.z());

    // evaluate gradient1 and gradient2 for katom:
    gradient1 = acos_prime *
                (-rij / (distij * distik) +
                 rij.dot(rik) * rik / (distij * distik * distik * distik)) *
                expij * expik;
    // gradient2
    gradient2 = ((rik / distik) * (-1.0 * gamma_sigma / (denomik * denomik))) *
                expij * expik;

    // add katom
    SP.AddToFitMatrix(A_, var,
                      least_sq_offset_ + 3 * nbeads_ * frame_counter_ + katom,
                      mpos, -gradient1.x(), -gradient2.x());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + nbeads_ + katom, mpos,
        -gradient1.y(), -gradient2.y());
    SP.AddToFitMatrix(
        A_, var,
        least_sq_offset_ + 3 * nbeads_ * frame_counter_ + 2 * nbeads_ + katom,
        mpos, -gradient1.z(), -gradient2.z());
  }
}
