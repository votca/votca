/* 
 * File:   main.cpp
 * Author: mashaya1
 *
 * Created on October 13, 2011, 10:52 PM
 */

#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <votca/csg/nblistgrid.h>
#include "csg_reupdate.h"
#include <votca/tools/linalg.h>

/*
 * 
 */
int main(int argc, char** argv) {

    CsgREupdate app;

    return app.Exec(argc, argv);
}

// program options are added here
void CsgREupdate::Initialize() {
    //cout << "Initializing csg_reupdate" << endl;
    
     CsgApplication::Initialize();
    // add RE options
    AddProgramOptions("RE Specific options")
    ("options", boost::program_options::value<string>(), 
            "  options file for coarse graining")
    ("gentable", boost::program_options::value<bool>(&_gentable)->default_value(false),
            "  only generate potential tables from given parameters, "
            "  NO RE update!");
    
    AddProgramOptions()
      ("top", boost::program_options::value<string > (), "  atomistic topology file (only needed for RE update)");
}

bool CsgREupdate::EvaluateOptions() {
    //cout << "Evaluating options" << endl;
    
    CsgApplication::EvaluateOptions();
    CheckRequired("options", "need to specify options file");
    if(!_gentable) {
        CheckRequired("trj", "no trajectory file specified");
        CheckRequired("top", "no topology file specified");
    }
    LoadOptions(OptionsMap()["options"].as<string>());

    return true;
}

void CsgREupdate::BeginEvaluate(Topology *top, Topology *top_atom){

    //cout << "Beginning evaluate" << endl;
    
    // initializing non-bonded interactions
    _nlamda = 0;
     for (list<Property*>::iterator iter = _nonbonded.begin();
             iter != _nonbonded.end(); ++iter) {

         string name = (*iter)->get("name").value();
         
         PotentialInfo *i = new PotentialInfo(_potentials.size(),
                                              false,
                                              _nlamda, *iter);

        // generate the bead lists
        BeadList beads1, beads2;

        beads1.Generate(*top, (*iter)->get("type1").value());
        beads2.Generate(*top, (*iter)->get("type2").value());

        if (beads1.size() == 0)
            throw std::runtime_error("Topology does not have beads of type \"" 
                    + (*iter)->get("type1").value() + "\"\n"
                "This was specified in type1 of interaction \"" + name + "\"");
        if (beads2.size() == 0)
            throw std::runtime_error("Topology does not have beads of type \"" 
                    + (*iter)->get("type2").value() + "\"\n"
                "This was specified in type2 of interaction \"" + name + "\"");
        // calculate normalization factor for rdf

        if ((*iter)->get("type1").value() == (*iter)->get("type2").value())
            i->rdf_norm = (beads1.size()*(beads2.size()) / 2.) / 
                    top->BoxVolume();
        else
            i->rdf_norm = (beads1.size() * beads2.size()) / top->BoxVolume();

        // let user know the properties of CG potential he/she selected
        cout << "We have " << i->potentialName << " CG potential" << endl;
        cout << "\t \t Between beads " << i->type1 << "-" << i->type2 << endl;
        cout << "\t \t With Function form " << i->potentialFunction << endl;
        cout << "\t \t And " << i->ucg->getOptParamSize() << " parameters to "
                "optimize" << endl;
        cout << "Potential range:" << endl;
        cout << "\t \t rmin    = " << i->rmin << " [nm]" << endl;
        cout << "\t \t rcutoff = " << i->rcut << " [nm]" << endl;
        cout << "\t \t step    = " << i->step   << " [nm]" << endl;

        // update parameter counter
        _nlamda += i->ucg->getOptParamSize();
        // add potential to container
        _potentials.push_back(i);
         
     }

     cout << "Total number of parameters to optimize: " << _nlamda << endl;

     _lamda.resize(_nlamda,false);
     // need to store initial guess of parameters in _lamda
    PotentialContainer::iterator potiter;

    for (potiter = _potentials.begin();
            potiter != _potentials.end(); ++potiter){

        int pos_start = (*potiter)->vec_pos;
        int pos_max  = pos_start + (*potiter)->ucg->getOptParamSize();

        for( int row = pos_start; row < pos_max; row++){

            int lamda_i = row - pos_start;
            _lamda(row) = (*potiter)->ucg->getOptParam(lamda_i);

        } // end row loop

    }// end potiter loop

    _dlamda.resize(_nlamda,false);
    _dlamda.clear();
    _DS.resize(_nlamda,false);
    _DS.clear();
    _HS.resize(_nlamda,_nlamda,false);
    _HS.clear();
    _dUFrame.resize(_nlamda,false);
    _dUFrame.clear();

    _nframes = 0.0; // no frames processed yet!
     // set Temperature
    _beta = _options.get("cg.inverse.kBT").as<double>();
    // relaxation parameter for update
    _relax = _options.get("cg.inverse.re.relax").as<double>();
     
    _UavgAA = 0.0;
    _UavgCG = 0.0;
    
}

void CsgREupdate::Run(){
    
    if( !_gentable )
    
        CsgApplication::Run();
    
    else { // only write potential tables for given parameters

        _nlamda = 0;
        for (list<Property*>::iterator iter = _nonbonded.begin();
                iter != _nonbonded.end(); ++iter) {

            string name = (*iter)->get("name").value();

            PotentialInfo *i = new PotentialInfo(_potentials.size(),
                    false,
                    _nlamda, *iter);

            // update parameter counter
            _nlamda += i->ucg->getOptParamSize();
            // add potential to container
            _potentials.push_back(i);

        }

        WriteOutFiles();

    }
    

}
void CsgREupdate::EndEvaluate(){
    
    //cout << "Ending RE Update" << endl;

    if (_nframes == 0 ){
        throw std::runtime_error("No frames to process! Please check your input.");
    }
    /* formulate _HS dlamda = - _DS   
     */

    REFormulateLinEq();

    REUpdateLamda();

    cout << "AA Ensemble Avg Energy :: " << _UavgAA << endl;
    cout << "CG Ensemble Avg Energy :: " << _UavgCG << endl;


    WriteOutFiles();
    
    cout <<"Finished RE update!\n";
}

void CsgREupdate::WriteOutFiles() {

    cout<< "Writing CG parameters and potential(s)\n";

    string potfile_extension = ".pot.new";
    string paramfile_extension = ".param.new";
    string file_name;

    
    
    PotentialContainer::iterator potiter;

    for (potiter = _potentials.begin();
            potiter != _potentials.end(); ++potiter) {

        //write potential table
        // construct meaningful outfile name
        file_name = (*potiter)->potentialName;
        file_name = file_name + potfile_extension;
        cout << "Writing file: " << file_name << endl;
        (*potiter)->ucg->SavePotTab(file_name,(*potiter)->step);

        //write parameter table
        // construct meaningful outfile name
        file_name = (*potiter)->potentialName;
        file_name = file_name + paramfile_extension;
        cout << "Writing file: " << file_name << endl;
        (*potiter)->ucg->SaveParam(file_name);
        
    }
    
}

// load user provided .xml option file
void CsgREupdate::LoadOptions(const string &file) {
    //cout << "Loading options" << endl;

    load_property_from_xml(_options, file);
    _nonbonded = _options.Select("cg.non-bonded");
    
}

// formulate _HS x = -_DS
void CsgREupdate::REFormulateLinEq() {

    /* compute CG ensemble avges of dU,d2U by dividing its
     * sum over all frames by total no. of frames
     */
    _DS /= ( (double)_nframes );
    _HS /= ( (double)_nframes );
    _UavgCG /= ( (double)_nframes );

    
    /* adding 4th term in eq. 52 of ref J. Chem. Phys. 134, 094112, 2011
     * to _HS
     */
    for( int row = 0; row < _nlamda; row++) {

        for( int col = 0; col < _nlamda; col++){

            _HS(row,col) += (-1.0*_DS(row)*_DS(col));
            // since at this step _DS(i) = -beta*<dU/dlamda_i>
        }// end loop over col
    } // end loop over row

    /* adding 1st term (i.e. aa ensemble avg) of eq. 51 to _DS
     * and of eq. 52 to _DH
     */
     PotentialContainer::iterator potiter;

     for (potiter = _potentials.begin();
             potiter != _potentials.end(); ++potiter){

         PotentialInfo *potinfo = *potiter;
         if( potinfo->bonded ){
             AAavgBonded(potinfo);
         } else {
             AAavgNonbonded(potinfo);
         }

     }
     
}

// update lamda = lamda + relax*dlamda
void CsgREupdate::REUpdateLamda() {

    // first solve _HS dx = -_DS
    ub::vector<double> residual(_nlamda);
    ub::vector<double> minusDS(_nlamda);

    minusDS = -_DS;

    try {
        votca::tools::linalg_cholesky_solve(_dlamda, _HS, minusDS);
    }
    catch (std::runtime_error){
        /* then can not use Newton-Raphson
         * go along steepest descent i.e. _dlamda = -_DS
         */
        ofstream out;
        string filename = "notsympos";
        out.open(filename.c_str());

        if(!out)
            throw runtime_error(string("error, cannot open file ") + filename);

        out << "**** Hessian NOT a positive definite! ****" << endl;
        out << "**** Taking steepest descent ****!" << endl;

        out.close();
        
        _dlamda = minusDS;
        
    }

    _lamda = _lamda + _relax * _dlamda ;

    // now update parameters of individual cg potentials
    PotentialContainer::iterator potiter;

    for (potiter = _potentials.begin();
            potiter != _potentials.end(); ++potiter){

        int pos_start = (*potiter)->vec_pos;
        int pos_max  = pos_start + (*potiter)->ucg->getOptParamSize();

        for( int row = pos_start; row < pos_max; row++){

            int lamda_i = row - pos_start;
            (*potiter)->ucg->setOptParam(lamda_i,_lamda(row));
            
        } // end row loop

    }// end potiter loop
    
}

// do non bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgNonbonded(PotentialInfo* potinfo) {

    int pos_start = potinfo->vec_pos;
    int pos_max   = pos_start + potinfo->ucg->getOptParamSize();
    int lamda_i,lamda_j;
    double dU_i,d2U_ij;
    double U;
    // compute avg AA energy
    U = 0.0;
    
    // assuming rdf bins are of same size
    double step = potinfo->aardf.x(2) - potinfo->aardf.x(1);
    
    for(int bin = 0; bin < potinfo->aardf.size(); bin++) {

        double r_hist = potinfo->aardf.x(bin);
        double r1 = r_hist - 0.5 * step;
        double r2 = r1 + step;
        double n_hist = potinfo->aardf.y(bin) * potinfo->rdf_norm *
                (4. / 3. * M_PI * (r2 * r2 * r2 - r1 * r1 * r1));
        U +=  n_hist *  potinfo->ucg->CalculateF(r_hist);

    }
    _UavgAA += U;
    
    // computing dU/dlamda and d2U/dlamda_i dlamda_j
    for( int row = pos_start; row < pos_max; row++){

        // ith parameter of this potential
        lamda_i = row - pos_start;

        // compute dU/dlamda and add to _DS
        dU_i = 0.0;
        for(int bin = 0; bin < potinfo->aardf.size(); bin++) {

            double r_hist = potinfo->aardf.x(bin);
            double r1 = r_hist - 0.5 * step;
            double r2 = r1 + step;
            double n_hist = potinfo->aardf.y(bin) * potinfo->rdf_norm *
                    (4. / 3. * M_PI * (r2 * r2 * r2 - r1 * r1 * r1));
            dU_i += n_hist * potinfo->ucg->CalculateDF(lamda_i,r_hist);
            
        } // end loop over hist
        _DS(row) += _beta * dU_i;
        
        for( int col = pos_start; col < pos_max; col++){

            lamda_j = col - pos_start;
            
            // compute d2U/dlamda_i dlamda_j and add to _HS
            d2U_ij = 0.0;
            for(int bin = 0; bin < potinfo->aardf.size(); bin++) {

                double r_hist = potinfo->aardf.x(bin);
                double r1 = r_hist - 0.5 * step;
                double r2 = r1 + step;
                double n_hist = potinfo->aardf.y(bin) * potinfo->rdf_norm *
                        (4. / 3. * M_PI * (r2 * r2 * r2 - r1 * r1 * r1));
                d2U_ij += n_hist *
                        potinfo->ucg->CalculateD2F(lamda_i, lamda_j, r_hist);

            } // end loop pair_iter

            _HS(row,col) += _beta * d2U_ij;
            
        } // end loop col

    } // end loop row
    
}

// do bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgBonded(PotentialInfo* potinfo) {

}

CsgApplication::Worker * CsgREupdate::ForkWorker(){

    CsgREupdateWorker *worker = new CsgREupdateWorker();

    // initialize worker
    worker->_options = _options;
    worker->_nonbonded = _nonbonded;
    worker->_nlamda = 0;

    for (list<Property*>::iterator iter = _nonbonded.begin();
             iter != _nonbonded.end(); ++iter) {

         PotentialInfo *i = new PotentialInfo(worker->_potentials.size(),
                                              false,
                                              worker->_nlamda, *iter);
        
        // update parameter counter
        worker->_nlamda += i->ucg->getOptParamSize();
        // add potential to container
        worker->_potentials.push_back(i);

     }

    worker->_lamda.resize(worker->_nlamda,false);
     // need to store initial guess of parameters in _lamda
    PotentialContainer::iterator potiter;

    for (potiter = worker->_potentials.begin();
            potiter != worker->_potentials.end(); ++potiter){

        int pos_start = (*potiter)->vec_pos;
        int pos_max  = pos_start + (*potiter)->ucg->getOptParamSize();

        for( int row = pos_start; row < pos_max; row++){

            int lamda_i = row - pos_start;
            worker->_lamda(row) = (*potiter)->ucg->getOptParam(lamda_i);

        } // end row loop

    }// end potiter loop

    worker->_DS.resize(worker->_nlamda,false);
    worker->_DS.clear();
    worker->_HS.resize(worker->_nlamda,worker->_nlamda,false);
    worker->_HS.clear();
    worker->_dUFrame.resize(worker->_nlamda,false);
    worker->_dUFrame.clear();

    worker->_nframes = 0.0; // no frames processed yet!
     // set Temperature
    worker->_beta = worker->_options.get("cg.inverse.kBT").as<double>();

    worker->_UavgAA = 0.0;
    worker->_UavgCG = 0.0;

    return worker;
    
}

void CsgREupdate::MergeWorker(Worker* worker) {

    CsgREupdateWorker *myCsgREupdateWorker;

    myCsgREupdateWorker = dynamic_cast<CsgREupdateWorker*> (worker);

    _UavgAA += myCsgREupdateWorker->_UavgAA;
    _UavgCG += myCsgREupdateWorker->_UavgCG;
    _nframes += myCsgREupdateWorker->_nframes;

    for(int row=0 ; row < _nlamda ; row++){
    
        _DS(row) += myCsgREupdateWorker->_DS(row);
        
        for( int col = 0; col < _nlamda; col++){
            _HS(row,col) += myCsgREupdateWorker->_HS(row,col);
        }
    }
    
}

void CsgREupdateWorker::EvalConfiguration(Topology *conf, Topology *conf_atom){

    /* as per Relative Entropy Ref.  J. Chem. Phys. 134, 094112, 2011
     * for each CG we need to compute dU/dlamda and d2U/dlamda_i dlamda_j
     * here, we add running sum of the second term of eq. 51 and store it
     * in _DS, ensemble averaging and first term addition is done in EndEvaluate
     * for eq. 52, here the running sum of second and third terms is computed
     * and added to _HS, addition of 1st and 4th term and ensemble average
     * is performed in EndEvalute
     */

    /* 3rd term of eq. 52 can not be computed in EvalBonded/EvalNonbonded
     * since only one CG potential is accessible in there.
     * hence store current frame dU/dlamda in _dUFrame!
     */

    // make sure _dUFrame is clear
    _dUFrame.clear();

    // loop over all CG potentials
     PotentialContainer::iterator potiter;

     for (potiter = _potentials.begin();
             potiter != _potentials.end(); ++potiter){

         PotentialInfo *potinfo = *potiter;
         if( potinfo->bonded ){
             EvalBonded(conf, potinfo);
         } else {
             EvalNonbonded(conf,potinfo);
         }

    }
    // update _DS and _HS
    for (int row = 0; row < _nlamda; row++) {

        _DS(row) += (-1.0 * _beta * _dUFrame(row));

        for (int col = 0; col < _nlamda; col++) {

            _HS(row, col) += (_beta * _beta * _dUFrame(row) * _dUFrame(col));

        }
    }

     _nframes++;
}


//do nonbonded potential related update stuff for the current frame in evalconfig
void CsgREupdateWorker::EvalNonbonded(Topology* conf, PotentialInfo* potinfo) {

    // generate the bead lists
    BeadList beads1, beads2;
    beads1.Generate(*conf, potinfo->type1);
    beads2.Generate(*conf, potinfo->type2);

    // check beads exists
    if(beads1.size() == 0)
            throw std::runtime_error("Topology does not have beads of type \""
                        + potinfo->type1 + "\"\n"
                    "This was specified in type1 of interaction \""
                    + potinfo->potentialName + "\"");
    if(beads2.size() == 0)
            throw std::runtime_error("Topology does not have beads of type \""
                        + potinfo->type2 + "\"\n"
                    "This was specified in type2 of interaction \""
                    + potinfo->potentialName + "\"");

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

    nb->setCutoff(potinfo->ucg->getCutOff());

    if( potinfo->type1 == potinfo->type2 ) { // same beads
        nb->Generate(beads1,true);
    } else { // different beads
        nb->Generate(beads1,beads2,true);
    }
    NBList::iterator pair_iter;

    int pos_start = potinfo->vec_pos;
    int pos_max = pos_start + potinfo->ucg->getOptParamSize();
    int lamda_i, lamda_j;
    double dU_i, d2U_ij;
    double U;

    // compute total energy
    U = 0.0;
    for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {
        U += potinfo->ucg->CalculateF((*pair_iter)->dist());
    }
    _UavgCG += U;

    // computing dU/dlamda and d2U/dlamda_i dlamda_j
    for (int row = pos_start; row < pos_max; row++) {

        // ith parameter of this potential
        lamda_i = row - pos_start;

        // compute dU/dlamda
        // double sum by iterating over all pairs
        dU_i = 0.0;
        for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {

            dU_i += potinfo->ucg->CalculateDF(lamda_i, (*pair_iter)->dist());

        } // end loop pair_iter

        _dUFrame(row) = dU_i;

        for (int col = pos_start; col < pos_max; col++) {

            lamda_j = col - pos_start;

            // compute d2U/dlamda_i dlamda_j and add to _HS
            // double sum by iterating over all pairs
            d2U_ij = 0.0;
            for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {

                d2U_ij += potinfo->ucg->CalculateD2F(lamda_i, lamda_j, (*pair_iter)->dist());

            } // end loop pair_iter

            _HS(row, col) += (-1.0 * _beta * d2U_ij);

        } // end loop col

    } // end loop row

    delete nb;
}

//do bonded potential related update stuff for the current frame in evalconfig
void CsgREupdateWorker::EvalBonded(Topology* conf, PotentialInfo* potinfo){

    // coming soon!

}

PotentialInfo::PotentialInfo(int index,
                             bool bonded_,
                             int vec_pos_, Property* options) {
    
    potentialIndex = index;
    _options = options;
    bonded   = bonded_;
    vec_pos  = vec_pos_;
    potentialName = options->get("name").value();
    type1 = options->get("type1").value();
    type2 = options->get("type2").value();
    potentialFunction = options->get("function").value();
    
    // compute output table grid points
    rmin = options->get("min").as<double>();
    rcut = options->get("max").as<double>();
    step = options->get("step").as<double>();
    
    // assign the user selected function form for this potential
    if( potentialFunction == "LJ126")
        ucg = new PotentialFunctionLJ126(rmin, rcut);
    else if (potentialFunction == "LJG")
        ucg = new PotentialFunctionLJG(rmin, rcut);
    else if (potentialFunction == "CBSPL"){

        // get number of B-splines coefficients which are to be optimized
        int nlam = options->get("nknots").as<int>();
        // get number of B-splines coefficients near cut-off which must be
        // fixed to zero to ensure potential and forces smoothly go to zero
        int ncut = options->get("nknots_cutoff").as<int>();

        ucg = new PotentialFunctionCBSPL(nlam, ncut, rmin, rcut);

    }
    else
        throw std::runtime_error("Function form \""
                        + potentialFunction + "\" selected for \""
                    + potentialName + "\" is not available yet.\n"
                    + "Please specify either \"LJ126, LJG, or CBSPL\" "
                    + "in options file.");

    // initialize cg potential with old parameters

    string oldparam_file_extension = ".param.cur";
    string oldparam_file_name = potentialName + oldparam_file_extension;

    ucg->setParam(oldparam_file_name);

    /* read/load reference AA CG-CG distribution
     */
    string aardf_file_extension = ".aa.rdf";
    string aardf_file_name = potentialName + aardf_file_extension;

    aardf.Load(aardf_file_name);

}
