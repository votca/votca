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
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/linalg.h>

// Local VOTCA includes
#include "votca/csg/nblistgrid.h"

// Local private VOTCA includes
#include "csg_reupdate.h"

using namespace std;

int main(int argc, char **argv) {
  CsgREupdate app;
  return app.Exec(argc, argv);
}

// program options are added here
void CsgREupdate::Initialize() {

  CsgApplication::Initialize();

  AddProgramOptions("RE Specific options")(
      "options", boost::program_options::value<string>(),
      "  options file for coarse graining")(
      "gentable",
      boost::program_options::value<bool>(&gentable_)->default_value(false),
      "  only generate potential tables from given parameters, "
      "  NO RE update!")("interaction", boost::program_options::value<string>(),
                         " [OPTIONAL] generate potential tables only for the "
                         "specified interactions, \n"
                         " only valid when 'gentable' is true")(
      "param-in-ext",
      boost::program_options::value<string>(&param_in_ext_)
          ->default_value("param.cur"),
      "  Extension of the input parameter tables")(
      "param-out-ext",
      boost::program_options::value<string>(&param_out_ext_)
          ->default_value("param.new"),
      "  Extension of the output parameter tables")(
      "pot-out-ext",
      boost::program_options::value<string>(&pot_out_ext_)
          ->default_value("pot.new"),
      "  Extension of the output potential tables")(
      "hessian-check",
      boost::program_options::value<bool>(&hessian_check_)->default_value(true),
      "  Disable the hessian check (mostly for testing)");
  AddProgramOptions()("top", boost::program_options::value<string>(),
                      "  atomistic topology file (only needed for RE update)");
}

bool CsgREupdate::EvaluateOptions() {

  CsgApplication::EvaluateOptions();
  CheckRequired("options", "need to specify options file");
  if (!gentable_) {
    CheckRequired("trj", "no trajectory file specified");
    CheckRequired("top", "no topology file specified");
  }
  LoadOptions(OptionsMap()["options"].as<string>());
  return true;
}

// load user provided .xml option file
void CsgREupdate::LoadOptions(const string &file) {

  options_.LoadFromXML(file);
  nonbonded_ = options_.Select("cg.non-bonded");
}

void CsgREupdate::BeginEvaluate(Topology *top, Topology *) {

  // initializing non-bonded interactions
  nlamda_ = 0;
  for (Property *prop : nonbonded_) {

    string name = prop->get("name").value();
    votca::Index id = potentials_.size();

    PotentialInfo *i =
        new PotentialInfo(id, false, nlamda_, param_in_ext_, prop, gentable_);

    Table *aardf = new Table();
    aardf->Load(name + ".dist.tgt");
    aardfs_.push_back(aardf);

    // generate the bead lists
    BeadList beads1, beads2;
    beads1.Generate(*top, prop->get("type1").value());
    beads2.Generate(*top, prop->get("type2").value());

    if (beads1.size() == 0) {
      throw std::runtime_error("Topology does not have beads of type \"" +
                               prop->get("type1").value() +
                               "\"\n"
                               "This was specified in type1 of interaction \"" +
                               name + "\"");
    }
    if (beads2.size() == 0) {
      throw std::runtime_error("Topology does not have beads of type \"" +
                               prop->get("type2").value() +
                               "\"\n"
                               "This was specified in type2 of interaction \"" +
                               name + "\"");
    }

    // calculate normalization factor for rdf
    double *rdfnorm = new double();

    if (prop->get("type1").value() == prop->get("type2").value()) {
      *rdfnorm =
          (double)(beads1.size() * beads2.size()) / 2. / top->BoxVolume();
    } else {
      *rdfnorm = (double)(beads1.size() * beads2.size()) / top->BoxVolume();
    }

    aardfnorms_.push_back(rdfnorm);

    // let user know the properties of CG potential he/she selected
    cout << "We have " << i->potentialName << " CG potential" << endl;
    cout << "\t \t Between beads " << i->type1 << "-" << i->type2 << endl;
    cout << "\t \t With Function form " << i->potentialFunction << endl;
    cout << "\t \t And " << i->ucg->getOptParamSize()
         << " parameters to "
            "optimize"
         << endl;
    cout << "Potential range:" << endl;
    cout << "\t \t rmin    = " << i->rmin << " [nm]" << endl;
    cout << "\t \t rcutoff = " << i->rcut << " [nm]" << endl;

    // update parameter counter
    nlamda_ += i->ucg->getOptParamSize();

    potentials_.push_back(i);
  }

  cout << "Total number of parameters to optimize: " << nlamda_ << endl;

  lamda_.resize(nlamda_);

  // need to store initial guess of parameters in  lamda_
  for (auto &potential_ : potentials_) {

    votca::Index pos_start = potential_->vec_pos;
    votca::Index pos_max = pos_start + potential_->ucg->getOptParamSize();

    for (votca::Index row = pos_start; row < pos_max; row++) {

      votca::Index lamda_i = row - pos_start;
      lamda_(row) = potential_->ucg->getOptParam(lamda_i);

    }  // end row loop

  }  // end potiter loop

  DS_ = Eigen::VectorXd::Zero(nlamda_);
  HS_ = Eigen::MatrixXd::Zero(nlamda_, nlamda_);
  dUFrame_ = Eigen::VectorXd::Zero(nlamda_);
  nframes_ = 0.0;  // no frames processed yet!

  // set Temperature
  beta_ = (1.0 / options_.get("cg.inverse.kBT").as<double>());

  // relaxation parameter for update
  relax_ = options_.get("cg.inverse.scale").as<double>();

  UavgAA_ = 0.0;
  UavgCG_ = 0.0;
}

void CsgREupdate::Run() {

  if (!gentable_) {
    CsgApplication::Run();
  } else {

    // only write potential tables for given parameters
    nlamda_ = 0;
    for (Property *prop : nonbonded_) {

      string name = prop->get("name").value();

      if (OptionsMap().count("interaction")) {

        Tokenizer tok(OptionsMap()["interaction"].as<string>(), ";");
        vector<string> vtok = tok.ToVector();
        vector<string>::iterator vtok_iter =
            find(vtok.begin(), vtok.end(), name);
        if (vtok_iter == vtok.end()) {
          continue;
        }
      }

      PotentialInfo *i = new PotentialInfo(potentials_.size(), false, nlamda_,
                                           param_in_ext_, prop, gentable_);

      // update parameter counter
      nlamda_ += i->ucg->getOptParamSize();

      potentials_.push_back(i);
    }

    WriteOutFiles();
  }
}

void CsgREupdate::EndEvaluate() {

  if (nframes_ == 0) {
    throw std::runtime_error("No frames to process! Please check your input.");
  }

  // formulate  HS_ dlamda = -  DS_

  REFormulateLinEq();

  cout << "Updating parameters" << endl;
  REUpdateLamda();

  cout << "AA Ensemble Avg Energy :: " << UavgAA_ << endl;
  cout << "CG Ensemble Avg Energy :: " << UavgCG_ << endl;

  WriteOutFiles();

  cout << "Finished RE update!\n";
}

void CsgREupdate::WriteOutFiles() {
  cout << "Writing CG parameters and potential(s)\n";
  string file_name;

  for (auto &potential_ : potentials_) {
    file_name = potential_->potentialName;
    file_name = file_name + "." + pot_out_ext_;
    cout << "Writing file: " << file_name << endl;
    potential_->ucg->SavePotTab(file_name,
                                potential_->options_->get("step").as<double>(),
                                potential_->options_->get("min").as<double>(),
                                potential_->options_->get("max").as<double>());
    // for gentable with no RE update no need to write-out parameters
    if (!gentable_) {
      file_name = potential_->potentialName;
      file_name = file_name + "." + param_out_ext_;
      cout << "Writing file: " << file_name << endl;
      potential_->ucg->SaveParam(file_name);
    }
  }
}

// formulate  HS_ x = - DS_
void CsgREupdate::REFormulateLinEq() {

  /* compute CG ensemble avges of dU,d2U by dividing its
   * sum over all frames by total no. of frames
   */

  DS_ /= ((double)nframes_);
  HS_ /= ((double)nframes_);
  UavgCG_ /= ((double)nframes_);

  /* adding 4th term in eq. 52 of ref J. Chem. Phys. 134, 094112, 2011
   * to  HS_
   */

  HS_ -= DS_ * DS_.transpose();

  /* adding 1st term (i.e. aa ensemble avg) of eq. 51 to  DS_
   * and of eq. 52 to  DH_
   */
  for (auto potinfo : potentials_) {

    if (potinfo->bonded) {
      AAavgBonded(potinfo);
    } else {
      AAavgNonbonded(potinfo);
    }
  }
}

// update lamda = lamda + relax*dlamda
void CsgREupdate::REUpdateLamda() {

  // first solve  HS_ dx = - DS_
  Eigen::LLT<Eigen::MatrixXd> cholesky(HS_);  // compute the Cholesky
                                              // decomposition of  HS_
  Eigen::VectorXd dlamda = cholesky.solve(-DS_);
  if (cholesky.info() == Eigen::ComputationInfo::NumericalIssue) {

    /* then can not use Newton-Raphson
     * steepest descent with line-search may be helpful
     * not implemented yet!
     */
    if (hessian_check_) {
      throw runtime_error(
          "Hessian NOT a positive definite!\n"
          "This can be a result of poor initial guess or "
          "ill-suited CG potential settings or poor CG sampling.\n");
    } else {

      cout << "****************** WARNING *****************" << endl;
      cout << "Hessian H NOT a positive definite!" << endl;
      cout << "This can be a result of poor initial guess or ill-suited CG "
              "potential settings or poor CG sampling."
           << endl;
      cout << "********************************************" << endl;
      cout << endl;
      cout << "However, you want to proceed despite non-positive definite H."
           << endl;
      cout << "Non-positive H does not gurantee relative entropy iterations to "
              "converge to minimum."
           << endl;
      cout << "You can turn on Hessian check by setting command line option "
              "--hessian-check=true."
           << endl;

      cout << "In this case, alternative update option is steepest descent."
           << endl;
      dlamda = -DS_;
    }
  }

  lamda_ = lamda_ + relax_ * dlamda;

  // now update parameters of individual cg potentials
  for (auto &potential_ : potentials_) {

    votca::Index pos_start = potential_->vec_pos;
    votca::Index pos_max = pos_start + potential_->ucg->getOptParamSize();

    for (votca::Index row = pos_start; row < pos_max; row++) {

      votca::Index lamda_i = row - pos_start;
      potential_->ucg->setOptParam(lamda_i, lamda_(row));

    }  // end row loop

  }  // end potiter loop
}

// do non bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgNonbonded(PotentialInfo *potinfo) {

  votca::Index pos_start = potinfo->vec_pos;
  votca::Index pos_max = pos_start + potinfo->ucg->getOptParamSize();
  double dU_i, d2U_ij;
  double U;
  votca::Index indx = potinfo->potentialIndex;

  // compute avg AA energy
  U = 0.0;

  // assuming rdf bins are of same size
  double step = aardfs_[indx]->x(2) - aardfs_[indx]->x(1);

  for (votca::Index bin = 0; bin < aardfs_[indx]->size(); bin++) {

    double r_hist = aardfs_[indx]->x(bin);
    double r1 = r_hist - 0.5 * step;
    double r2 = r1 + step;
    double n_hist =
        aardfs_[indx]->y(bin) * (*aardfnorms_[indx]) *
        (4. / 3. * votca::tools::conv::Pi * (r2 * r2 * r2 - r1 * r1 * r1));

    if (n_hist > 0.0) {
      U += n_hist * potinfo->ucg->CalculateF(r_hist);
    }
  }

  UavgAA_ += U;

  // computing dU/dlamda and d2U/dlamda_i dlamda_j
  for (votca::Index row = pos_start; row < pos_max; row++) {

    // ith parameter of this potential
    votca::Index lamda_i = row - pos_start;

    // compute dU/dlamda and add to  DS_
    dU_i = 0.0;

    for (votca::Index bin = 0; bin < aardfs_[indx]->size(); bin++) {

      double r_hist = aardfs_[indx]->x(bin);
      double r1 = r_hist - 0.5 * step;
      double r2 = r1 + step;
      double n_hist =
          aardfs_[indx]->y(bin) * (*aardfnorms_[indx]) *
          (4. / 3. * votca::tools::conv::Pi * (r2 * r2 * r2 - r1 * r1 * r1));

      if (n_hist > 0.0) {
        dU_i += n_hist * potinfo->ucg->CalculateDF(lamda_i, r_hist);
      }

    }  // end loop over hist

    DS_(row) += (beta_ * dU_i);

    for (votca::Index col = row; col < pos_max; col++) {

      votca::Index lamda_j = col - pos_start;

      // compute d2U/dlamda_i dlamda_j and add to  HS_
      d2U_ij = 0.0;

      for (votca::Index bin = 0; bin < aardfs_[indx]->size(); bin++) {

        double r_hist = aardfs_[indx]->x(bin);
        double r1 = r_hist - 0.5 * step;
        double r2 = r1 + step;
        double n_hist =
            aardfs_[indx]->y(bin) * (*aardfnorms_[indx]) *
            (4. / 3. * votca::tools::conv::Pi * (r2 * r2 * r2 - r1 * r1 * r1));

        if (n_hist > 0.0) {
          d2U_ij +=
              n_hist * potinfo->ucg->CalculateD2F(lamda_i, lamda_j, r_hist);
        }

      }  // end loop pair_iter

      HS_(row, col) += (beta_ * d2U_ij);
      if (row != col) {
        HS_(col, row) += (beta_ * d2U_ij);
      }
    }  // end loop col

  }  // end loop row
}

// do bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgBonded(PotentialInfo *) {}

std::unique_ptr<CsgApplication::Worker> CsgREupdate::ForkWorker() {

  auto worker = std::make_unique<CsgREupdateWorker>();

  // initialize worker
  worker->options_ = options_;
  worker->nonbonded_ = nonbonded_;
  worker->nlamda_ = 0;

  for (Property *prop : nonbonded_) {

    PotentialInfo *i = new PotentialInfo(worker->potentials_.size(), false,
                                         worker->nlamda_, param_in_ext_, prop);
    // update parameter counter
    worker->nlamda_ += i->ucg->getOptParamSize();

    worker->potentials_.push_back(i);
  }

  worker->lamda_.resize(worker->nlamda_);

  // need to store initial guess of parameters in  lamda_
  for (PotentialInfo *pot : worker->potentials_) {

    votca::Index pos_start = pot->vec_pos;
    votca::Index pos_max = pos_start + pot->ucg->getOptParamSize();

    for (votca::Index row = pos_start; row < pos_max; row++) {

      votca::Index lamda_i = row - pos_start;
      worker->lamda_(row) = pot->ucg->getOptParam(lamda_i);

    }  // end row loop

  }  // end potiter loop

  worker->DS_ = Eigen::VectorXd::Zero(worker->nlamda_);
  worker->HS_ = Eigen::MatrixXd::Zero(worker->nlamda_, worker->nlamda_);
  worker->dUFrame_ = Eigen::VectorXd::Zero(worker->nlamda_);
  worker->nframes_ = 0.0;  // no frames processed yet!
  // set Temperature
  worker->beta_ = (1.0 / worker->options_.get("cg.inverse.kBT").as<double>());
  worker->UavgCG_ = 0.0;

  return worker;
}

void CsgREupdate::MergeWorker(Worker *worker) {

  CsgREupdateWorker *myCsgREupdateWorker;
  myCsgREupdateWorker = dynamic_cast<CsgREupdateWorker *>(worker);
  UavgCG_ += myCsgREupdateWorker->UavgCG_;
  nframes_ += myCsgREupdateWorker->nframes_;
  DS_ += myCsgREupdateWorker->DS_;
  HS_ += myCsgREupdateWorker->HS_;
}

void CsgREupdateWorker::EvalConfiguration(Topology *conf, Topology *) {

  /* as per Relative Entropy Ref.  J. Chem. Phys. 134, 094112, 2011
   * for each CG we need to compute dU/dlamda and d2U/dlamda_i dlamda_j
   * here, we add running sum of the second term of eq. 51 and store it
   * in  DS_, ensemble averaging and first term addition is done in EndEvaluate
   * for eq. 52, here the running sum of second and third terms is computed
   * and added to  HS_, addition of 1st and 4th term and ensemble average
   * is performed in EndEvalute
   * 3rd term of eq. 52 can not be computed in EvalBonded/EvalNonbonded
   * since only one CG potential is accessible in there.
   * hence store current frame dU/dlamda in  dUFrame_!
   */

  dUFrame_.setZero();

  for (PotentialInfo *potinfo : potentials_) {

    if (potinfo->bonded) {
      EvalBonded(conf, potinfo);
    } else {
      EvalNonbonded(conf, potinfo);
    }
  }
  // update  DS_ and  HS_
  DS_ -= beta_ * dUFrame_;
  HS_ += beta_ * beta_ * dUFrame_ * dUFrame_.transpose();

  nframes_++;
}

// do nonbonded potential related update stuff for the current frame in
// evalconfig
void CsgREupdateWorker::EvalNonbonded(Topology *conf, PotentialInfo *potinfo) {

  BeadList beads1, beads2;
  beads1.Generate(*conf, potinfo->type1);
  beads2.Generate(*conf, potinfo->type2);

  if (beads1.size() == 0) {
    throw std::runtime_error("Topology does not have beads of type \"" +
                             potinfo->type1 +
                             "\"\n"
                             "This was specified in type1 of interaction \"" +
                             potinfo->potentialName + "\"");
  }
  if (beads2.size() == 0) {
    throw std::runtime_error("Topology does not have beads of type \"" +
                             potinfo->type2 +
                             "\"\n"
                             "This was specified in type2 of interaction \"" +
                             potinfo->potentialName + "\"");
  }

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

  nb->setCutoff(potinfo->ucg->getCutOff());

  if (potinfo->type1 == potinfo->type2) {  // same beads
    nb->Generate(beads1, true);
  } else {  // different beads
    nb->Generate(beads1, beads2, true);
  }

  votca::Index pos_start = potinfo->vec_pos;
  votca::Index pos_max = pos_start + potinfo->ucg->getOptParamSize();

  // compute total energy
  double U = 0.0;
  for (auto &pair_iter : *nb) {
    U += potinfo->ucg->CalculateF(pair_iter->dist());
  }

  UavgCG_ += U;

  // computing dU/dlamda and d2U/dlamda_i dlamda_j
  for (votca::Index row = pos_start; row < pos_max; row++) {

    votca::Index lamda_i = row - pos_start;

    double dU_i = 0.0;
    for (auto &pair_iter : *nb) {
      dU_i += potinfo->ucg->CalculateDF(lamda_i, pair_iter->dist());
    }

    dUFrame_(row) = dU_i;

    for (votca::Index col = row; col < pos_max; col++) {

      votca::Index lamda_j = col - pos_start;
      double d2U_ij = 0.0;

      for (auto &pair_iter : *nb) {
        d2U_ij +=
            potinfo->ucg->CalculateD2F(lamda_i, lamda_j, pair_iter->dist());
      }

      HS_(row, col) -= beta_ * d2U_ij;
      if (row != col) {
        HS_(col, row) -= beta_ * d2U_ij;
      }
    }  // end loop col

  }  // end loop row
}

// do bonded potential related update stuff for the current frame in evalconfig
void CsgREupdateWorker::EvalBonded(Topology *, PotentialInfo *) {}

PotentialInfo::PotentialInfo(votca::Index index, bool bonded_,
                             votca::Index vec_pos_, string &param_in_ext_,
                             Property *options, bool gentable) {

  potentialIndex = index;
  bonded = bonded_;
  vec_pos = vec_pos_;
  options_ = options;

  potentialName = options_->get("name").value();
  type1 = options_->get("type1").value();
  type2 = options_->get("type2").value();
  potentialFunction = options_->get("re.function").value();

  rmin = options_->get("min").as<double>();
  rcut = options_->get("max").as<double>();

  // assign the user selected function form for this potential
  if (potentialFunction == "lj126") {
    ucg = new PotentialFunctionLJ126(potentialName, rmin, rcut);
  } else if (potentialFunction == "ljg") {
    ucg = new PotentialFunctionLJG(potentialName, rmin, rcut);
  } else if (potentialFunction == "cbspl") {
    // get number of B-splines coefficients which are to be optimized
    votca::Index nlam = options_->get("re.cbspl.nknots").as<votca::Index>();

    if (!gentable) {
      // determine minimum for B-spline from CG-MD rdf

      Table dist;
      string filename = potentialName + ".dist.new";
      try {
        dist.Load(filename);
        // for( Index i = 0; i < dist.size(); i++)
        // In the highly repulsive core (<0.1 nm) at some locations
        // an unphysical non-zero RDF value may occur,
        // so it would be better to estimate Rmin loop
        // through RDF values from Rcut to zero instead of zero to Rcut.
        double new_min = 0.0;
        for (votca::Index i = dist.size() - 2; i > 0; i--) {
          if (dist.y(i) < 1.0e-4) {
            new_min = dist.x(i + 1);
            break;
          }
        }
        // setting extrapolation region based on CG rdf
        rmin = new_min;

      } catch (std::runtime_error &) {
        throw std::runtime_error(
            "Missing file for CG rdf for the interaction " + potentialName +
            "." + " Make sure CG RDF for the pair " + potentialName +
            " is computed after CG-MD simulation.");
      }
    }
    ucg = new PotentialFunctionCBSPL(potentialName, nlam, rmin, rcut);
  } else {
    throw std::runtime_error(
        "Function form \"" + potentialFunction + "\" selected for \"" +
        potentialName + "\" is not available yet.\n" +
        "Please specify either \"lj126, ljg, or cbspl\" " + "in options file.");
  }
  // initialize cg potential with old parameters
  string oldparam_file_name = potentialName + "." + param_in_ext_;
  ucg->setParam(oldparam_file_name);
}
