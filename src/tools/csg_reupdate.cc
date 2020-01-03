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

#include "csg_reupdate.h"
#include "../../include/votca/csg/nblistgrid.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <votca/tools/constants.h>
#include <votca/tools/linalg.h>
/*
 *
 */

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
      boost::program_options::value<bool>(&_gentable)->default_value(false),
      "  only generate potential tables from given parameters, "
      "  NO RE update!")("interaction", boost::program_options::value<string>(),
                         " [OPTIONAL] generate potential tables only for the "
                         "specified interactions, \n"
                         " only valid when 'gentable' is true")(
      "param-in-ext",
      boost::program_options::value<string>(&_param_in_ext)
          ->default_value("param.cur"),
      "  Extension of the input parameter tables")(
      "param-out-ext",
      boost::program_options::value<string>(&_param_out_ext)
          ->default_value("param.new"),
      "  Extension of the output parameter tables")(
      "pot-out-ext",
      boost::program_options::value<string>(&_pot_out_ext)
          ->default_value("pot.new"),
      "  Extension of the output potential tables")(
      "hessian-check",
      boost::program_options::value<bool>(&_hessian_check)->default_value(true),
      "  Extension of the output potential tables");
  AddProgramOptions()("top", boost::program_options::value<string>(),
                      "  atomistic topology file (only needed for RE update)");
}

bool CsgREupdate::EvaluateOptions() {

  CsgApplication::EvaluateOptions();
  CheckRequired("options", "need to specify options file");
  if (!_gentable) {
    CheckRequired("trj", "no trajectory file specified");
    CheckRequired("top", "no topology file specified");
  }
  LoadOptions(OptionsMap()["options"].as<string>());
  return true;
}

// load user provided .xml option file
void CsgREupdate::LoadOptions(const string &file) {

  _options.LoadFromXML(file);
  _nonbonded = _options.Select("cg.non-bonded");
}

void CsgREupdate::BeginEvaluate(Topology *top, Topology *) {

  // initializing non-bonded interactions
  _nlamda = 0;
  for (Property *prop : _nonbonded) {

    string name = prop->get("name").value();
    votca::Index id = _potentials.size();

    PotentialInfo *i =
        new PotentialInfo(id, false, _nlamda, _param_in_ext, prop, _gentable);

    Table *aardf = new Table();
    aardf->Load(name + ".dist.tgt");
    _aardfs.push_back(aardf);

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

    _aardfnorms.push_back(rdfnorm);

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
    _nlamda += i->ucg->getOptParamSize();

    _potentials.push_back(i);
  }

  cout << "Total number of parameters to optimize: " << _nlamda << endl;

  _lamda.resize(_nlamda);

  // need to store initial guess of parameters in _lamda
  for (auto &_potential : _potentials) {

    votca::Index pos_start = _potential->vec_pos;
    votca::Index pos_max = pos_start + _potential->ucg->getOptParamSize();

    for (votca::Index row = pos_start; row < pos_max; row++) {

      votca::Index lamda_i = row - pos_start;
      _lamda(row) = _potential->ucg->getOptParam(lamda_i);

    }  // end row loop

  }  // end potiter loop

  _DS = Eigen::VectorXd::Zero(_nlamda);
  _HS = Eigen::MatrixXd::Zero(_nlamda, _nlamda);
  _dUFrame = Eigen::VectorXd::Zero(_nlamda);
  _nframes = 0.0;  // no frames processed yet!

  // set Temperature
  _beta = (1.0 / _options.get("cg.inverse.kBT").as<double>());

  // relaxation parameter for update
  _relax = _options.get("cg.inverse.scale").as<double>();

  // whether to take steepest descent in case of non-symmetric positive H
  //_dosteep = _options.get("cg.inverse.re.do_steep").as<bool>();

  _UavgAA = 0.0;
  _UavgCG = 0.0;
}

void CsgREupdate::Run() {

  if (!_gentable) {
    CsgApplication::Run();
  } else {

    // only write potential tables for given parameters
    _nlamda = 0;
    for (Property *prop : _nonbonded) {

      string name = prop->get("name").value();

      if (OptionsMap().count("interaction")) {

        Tokenizer tok(_op_vm["interaction"].as<string>(), ";");
        vector<string> vtok = tok.ToVector();
        vector<string>::iterator vtok_iter =
            find(vtok.begin(), vtok.end(), name);
        if (vtok_iter == vtok.end()) {
          continue;
        }
      }

      PotentialInfo *i = new PotentialInfo(_potentials.size(), false, _nlamda,
                                           _param_in_ext, prop, _gentable);

      // update parameter counter
      _nlamda += i->ucg->getOptParamSize();

      _potentials.push_back(i);
    }

    WriteOutFiles();
  }
}

void CsgREupdate::EndEvaluate() {

  if (_nframes == 0) {
    throw std::runtime_error("No frames to process! Please check your input.");
  }

  // formulate _HS dlamda = - _DS

  REFormulateLinEq();

  cout << "Updating parameters" << endl;
  REUpdateLamda();

  cout << "AA Ensemble Avg Energy :: " << _UavgAA << endl;
  cout << "CG Ensemble Avg Energy :: " << _UavgCG << endl;

  WriteOutFiles();

  cout << "Finished RE update!\n";
}

void CsgREupdate::WriteOutFiles() {
  cout << "Writing CG parameters and potential(s)\n";
  string file_name;

  for (auto &_potential : _potentials) {
    file_name = _potential->potentialName;
    file_name = file_name + "." + _pot_out_ext;
    cout << "Writing file: " << file_name << endl;
    _potential->ucg->SavePotTab(file_name,
                                _potential->_options->get("step").as<double>(),
                                _potential->_options->get("min").as<double>(),
                                _potential->_options->get("max").as<double>());
    // for gentable with no RE update no need to write-out parameters
    if (!_gentable) {
      file_name = _potential->potentialName;
      file_name = file_name + "." + _param_out_ext;
      cout << "Writing file: " << file_name << endl;
      _potential->ucg->SaveParam(file_name);
    }
  }
}

// formulate _HS x = -_DS
void CsgREupdate::REFormulateLinEq() {

  /* compute CG ensemble avges of dU,d2U by dividing its
   * sum over all frames by total no. of frames
   */

  _DS /= ((double)_nframes);
  _HS /= ((double)_nframes);
  _UavgCG /= ((double)_nframes);

  /* adding 4th term in eq. 52 of ref J. Chem. Phys. 134, 094112, 2011
   * to _HS
   */
  for (votca::Index row = 0; row < _nlamda; row++) {

    for (votca::Index col = row; col < _nlamda; col++) {

      _HS(row, col) += (-1.0 * _DS(row) * _DS(col));
      _HS(col, row) += (-1.0 * _DS(row) * _DS(col));
      // since at this step _DS(i) = -beta*<dU/dlamda_i>cg

    }  // end loop over col

  }  // end loop over row

  /* adding 1st term (i.e. aa ensemble avg) of eq. 51 to _DS
   * and of eq. 52 to _DH
   */
  for (auto potinfo : _potentials) {

    if (potinfo->bonded) {
      AAavgBonded(potinfo);
    } else {
      AAavgNonbonded(potinfo);
    }
  }
}

// update lamda = lamda + relax*dlamda
void CsgREupdate::REUpdateLamda() {

  // first solve _HS dx = -_DS
  Eigen::LLT<Eigen::MatrixXd> cholesky(_HS);  // compute the Cholesky
                                              // decomposition of _HS
  Eigen::VectorXd dlamda = cholesky.solve(-_DS);
  if (cholesky.info() == Eigen::ComputationInfo::NumericalIssue) {

    /* then can not use Newton-Raphson
     * steepest descent with line-search may be helpful
     * not implemented yet!
     */
    // if(!_dosteep){
    if (_hessian_check) {
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
      dlamda = -_DS;
    }
  }

  _lamda = _lamda + _relax * dlamda;

  // now update parameters of individual cg potentials
  for (auto &_potential : _potentials) {

    votca::Index pos_start = _potential->vec_pos;
    votca::Index pos_max = pos_start + _potential->ucg->getOptParamSize();

    for (votca::Index row = pos_start; row < pos_max; row++) {

      votca::Index lamda_i = row - pos_start;
      _potential->ucg->setOptParam(lamda_i, _lamda(row));

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
  double step = _aardfs[indx]->x(2) - _aardfs[indx]->x(1);

  for (votca::Index bin = 0; bin < _aardfs[indx]->size(); bin++) {

    double r_hist = _aardfs[indx]->x(bin);
    double r1 = r_hist - 0.5 * step;
    double r2 = r1 + step;
    double n_hist =
        _aardfs[indx]->y(bin) * (*_aardfnorms[indx]) *
        (4. / 3. * votca::tools::conv::Pi * (r2 * r2 * r2 - r1 * r1 * r1));

    if (n_hist > 0.0) {
      U += n_hist * potinfo->ucg->CalculateF(r_hist);
    }
  }

  _UavgAA += U;

  // computing dU/dlamda and d2U/dlamda_i dlamda_j
  for (votca::Index row = pos_start; row < pos_max; row++) {

    // ith parameter of this potential
    votca::Index lamda_i = row - pos_start;

    // compute dU/dlamda and add to _DS
    dU_i = 0.0;

    for (votca::Index bin = 0; bin < _aardfs[indx]->size(); bin++) {

      double r_hist = _aardfs[indx]->x(bin);
      double r1 = r_hist - 0.5 * step;
      double r2 = r1 + step;
      double n_hist =
          _aardfs[indx]->y(bin) * (*_aardfnorms[indx]) *
          (4. / 3. * votca::tools::conv::Pi * (r2 * r2 * r2 - r1 * r1 * r1));

      if (n_hist > 0.0) {
        dU_i += n_hist * potinfo->ucg->CalculateDF(lamda_i, r_hist);
      }

    }  // end loop over hist

    _DS(row) += (_beta * dU_i);

    for (votca::Index col = row; col < pos_max; col++) {

      votca::Index lamda_j = col - pos_start;

      // compute d2U/dlamda_i dlamda_j and add to _HS
      d2U_ij = 0.0;

      for (votca::Index bin = 0; bin < _aardfs[indx]->size(); bin++) {

        double r_hist = _aardfs[indx]->x(bin);
        double r1 = r_hist - 0.5 * step;
        double r2 = r1 + step;
        double n_hist =
            _aardfs[indx]->y(bin) * (*_aardfnorms[indx]) *
            (4. / 3. * votca::tools::conv::Pi * (r2 * r2 * r2 - r1 * r1 * r1));

        if (n_hist > 0.0) {
          d2U_ij +=
              n_hist * potinfo->ucg->CalculateD2F(lamda_i, lamda_j, r_hist);
        }

      }  // end loop pair_iter

      _HS(row, col) += (_beta * d2U_ij);
      _HS(col, row) += (_beta * d2U_ij);

    }  // end loop col

  }  // end loop row
}

// do bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgBonded(PotentialInfo *) {}

CsgApplication::Worker *CsgREupdate::ForkWorker() {

  CsgREupdateWorker *worker = new CsgREupdateWorker();

  // initialize worker
  worker->_options = _options;
  worker->_nonbonded = _nonbonded;
  worker->_nlamda = 0;

  for (Property *prop : _nonbonded) {

    PotentialInfo *i = new PotentialInfo(worker->_potentials.size(), false,
                                         worker->_nlamda, _param_in_ext, prop);
    // update parameter counter
    worker->_nlamda += i->ucg->getOptParamSize();

    worker->_potentials.push_back(i);
  }

  worker->_lamda.resize(worker->_nlamda);

  // need to store initial guess of parameters in _lamda
  for (PotentialInfo *pot : worker->_potentials) {

    votca::Index pos_start = pot->vec_pos;
    votca::Index pos_max = pos_start + pot->ucg->getOptParamSize();

    for (votca::Index row = pos_start; row < pos_max; row++) {

      votca::Index lamda_i = row - pos_start;
      worker->_lamda(row) = pot->ucg->getOptParam(lamda_i);

    }  // end row loop

  }  // end potiter loop

  worker->_DS = Eigen::VectorXd::Zero(worker->_nlamda);
  worker->_HS = Eigen::MatrixXd::Zero(worker->_nlamda, worker->_nlamda);
  worker->_dUFrame = Eigen::VectorXd::Zero(worker->_nlamda);
  worker->_nframes = 0.0;  // no frames processed yet!
  // set Temperature
  worker->_beta = (1.0 / worker->_options.get("cg.inverse.kBT").as<double>());
  worker->_UavgCG = 0.0;

  return worker;
}

void CsgREupdate::MergeWorker(Worker *worker) {

  CsgREupdateWorker *myCsgREupdateWorker;
  myCsgREupdateWorker = dynamic_cast<CsgREupdateWorker *>(worker);
  _UavgCG += myCsgREupdateWorker->_UavgCG;
  _nframes += myCsgREupdateWorker->_nframes;
  _DS += myCsgREupdateWorker->_DS;
  _HS += myCsgREupdateWorker->_HS;
}

void CsgREupdateWorker::EvalConfiguration(Topology *conf, Topology *) {

  /* as per Relative Entropy Ref.  J. Chem. Phys. 134, 094112, 2011
   * for each CG we need to compute dU/dlamda and d2U/dlamda_i dlamda_j
   * here, we add running sum of the second term of eq. 51 and store it
   * in _DS, ensemble averaging and first term addition is done in EndEvaluate
   * for eq. 52, here the running sum of second and third terms is computed
   * and added to _HS, addition of 1st and 4th term and ensemble average
   * is performed in EndEvalute
   * 3rd term of eq. 52 can not be computed in EvalBonded/EvalNonbonded
   * since only one CG potential is accessible in there.
   * hence store current frame dU/dlamda in _dUFrame!
   */

  _dUFrame.setZero();

  for (PotentialInfo *potinfo : _potentials) {

    if (potinfo->bonded) {
      EvalBonded(conf, potinfo);
    } else {
      EvalNonbonded(conf, potinfo);
    }
  }
  // update _DS and _HS
  _DS -= _beta * _dUFrame;
  _HS += _beta * _beta * _dUFrame * _dUFrame.transpose();

  _nframes++;
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

  if (_options.exists("cg.nbsearch")) {

    if (_options.get("cg.nbsearch").as<string>() == "grid") {
      gridsearch = true;
    } else if (_options.get("cg.nbsearch").as<string>() == "simple") {
      gridsearch = false;
    } else {
      throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
  }

  if (gridsearch) {
    nb = std::make_unique<NBList>(NBListGrid());
  } else {
    nb = std::make_unique<NBList>(NBList());
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

  _UavgCG += U;

  // computing dU/dlamda and d2U/dlamda_i dlamda_j
  for (votca::Index row = pos_start; row < pos_max; row++) {

    votca::Index lamda_i = row - pos_start;

    double dU_i = 0.0;
    for (auto &pair_iter : *nb) {
      dU_i += potinfo->ucg->CalculateDF(lamda_i, pair_iter->dist());
    }

    _dUFrame(row) = dU_i;

    for (votca::Index col = row; col < pos_max; col++) {

      votca::Index lamda_j = col - pos_start;
      double d2U_ij = 0.0;

      for (auto &pair_iter : *nb) {
        d2U_ij +=
            potinfo->ucg->CalculateD2F(lamda_i, lamda_j, pair_iter->dist());
      }

      _HS(row, col) -= _beta * d2U_ij;
      _HS(col, row) -= _beta * d2U_ij;
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
  _options = options;

  potentialName = _options->get("name").value();
  type1 = _options->get("type1").value();
  type2 = _options->get("type2").value();
  potentialFunction = _options->get("re.function").value();

  rmin = _options->get("min").as<double>();
  rcut = _options->get("max").as<double>();

  // assign the user selected function form for this potential
  if (potentialFunction == "lj126") {
    ucg = new PotentialFunctionLJ126(potentialName, rmin, rcut);
  } else if (potentialFunction == "ljg") {
    ucg = new PotentialFunctionLJG(potentialName, rmin, rcut);
  } else if (potentialFunction == "cbspl") {
    // get number of B-splines coefficients which are to be optimized
    votca::Index nlam = _options->get("re.cbspl.nknots").as<votca::Index>();

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
