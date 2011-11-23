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
#include <votca/tools/linalg.h>

#include "csg_reupdate.h"

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
    ("genref",boost::program_options::value<bool>(&_genref)->default_value(false),
            "generate referece All-atom neigborlist histogram");
    
}

bool CsgREupdate::EvaluateOptions() {
    //cout << "Evaluating options" << endl;
    
    CsgApplication::EvaluateOptions();
    CheckRequired("options", "need to specify options file");
    CheckRequired("trj", "no trajectory file specified");
    LoadOptions(OptionsMap()["options"].as<string>());

    return true;
}

void CsgREupdate::BeginEvaluate(Topology *top, Topology *top_atom){

    //cout << "Beginning evaluate" << endl;
    if(_genref)
        cout << "NO RE update! only process reference AA trajectories" << endl;
    
    // initializing non-bonded interactions
    _nlamda = 0;
     for (list<Property*>::iterator iter = _nonbonded.begin();
             iter != _nonbonded.end(); ++iter) {

         PotentialInfo *i = new PotentialInfo(top, _potentials.size(),
                                              false, _genref,
                                              _nlamda, *iter);
         // update parameter counter
         _nlamda += i->ucg->getParamSize();
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
        int pos_max  = pos_start + (*potiter)->ucg->getParamSize();

        for( int row = pos_start; row < pos_max; row++){

            int lamda_i = row - pos_start;
            _lamda(row) = (*potiter)->ucg->getParam(lamda_i);

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
    _relax = _options.get("cg.inverse.re.relax").as<double>();
     
    _UavgAA = 0.0;
    _UavgCG = 0.0;
    _HS_Scale = 0.0;
    
     if(_genref){ // then we need to write mapped cg topology for a single frame

         string out = "confout.gro";
         cout << "writing coarse-grained trajectory to " << out << endl;

         _writer = TrjWriterFactory().Create(out);
         if (_writer == NULL)
             throw runtime_error("output format not supported: " + out);

         _writer->Open(out);
         
     }
}

void CsgREupdate::EndEvaluate(){
    
    //cout << "Ending RE Update" << endl;

    if (_nframes == 0 ){
        throw std::runtime_error("No frames to process! Please check your input.");
    }
    /* formulate _HS dlamda = - _DS 
     * and solve it only if more than 0 frames processed     
     */
    if( !_genref ){
    
        REFormulateLinEq();

        REUpdateLamda();

        cout << "AA Ensemble Avg Energy :: " << _UavgAA << endl;
        cout << "CG Ensemble Avg Energy :: " << _UavgCG << endl;
        
    } else {
        // compute reference aa avg histograms
        AAavgHist();

        // close and delete cg map topology writer
         _writer->Close();
         delete _writer;
      
    }
    
    WriteOutFiles();
    
    cout <<"Finished RE update!\n";
}

void CsgREupdate::WriteOutFiles() {

    //cout << "Processed total of " << _nframes << " frames" << endl;
    cout<< "Writing CG parameters and potential(s)\n";

    string potfile_extension = ".pot.new";
    string paramfile_extension = ".param.new";
    string aahistfile_extension = ".aa.nbhist";
    string file_name;
    Table pot_tab;
    Table param_tab;
    Table aahist_tab;
    // table with error column or no error column
    pot_tab.SetHasYErr(false);
    param_tab.SetHasYErr(true);
    aahist_tab.SetHasYErr(false);
    
    PotentialContainer::iterator potiter;;

    for (potiter = _potentials.begin();
            potiter != _potentials.end(); ++potiter) {

        //write potential table
        // construct meaningful outfile name
        file_name = (*potiter)->potentialName;
        file_name = file_name + potfile_extension;

        // resize table
        int ngrid = (*potiter)->pottblgrid.size();
        pot_tab.resize(ngrid, false);

        // print output file names on stdout
        cout << "Writing file: " << file_name << endl;

        // loop over output grid points
        for (int i = 0; i < ngrid; i++) {

            double out_x = (*potiter)->pottblgrid(i);
            double out_u = (*potiter)->ucg->CalculateF(out_x);
            // put point, result, flag at point out_x into the table
            pot_tab.set(i, out_x, out_u , 'i');

        }
        // save table in the file
        pot_tab.Save(file_name);
        // clear the table for the next potential
        pot_tab.clear();

        //write parameter table
        // construct meaningful outfile name
        file_name = (*potiter)->potentialName;
        file_name = file_name + paramfile_extension;

        // resize table
        int nparam = (*potiter)->ucg->getParamSize();
        param_tab.resize(nparam, false);

        // print output file names on stdout
        cout << "Writing file: " << file_name << endl;

        // loop over output paramaeters
        for (int i = 0; i < nparam; i++) {

            int out_indx  = i;
            double out_param = (*potiter)->ucg->getParam(out_indx);
            int row = (*potiter)->vec_pos + i;
            // put point, result, flag and error at point out_x into the table
            param_tab.set(i, out_indx, out_param, 'i', _dlamda(row));

        }
        // save table in the file
        param_tab.Save(file_name);
        // clear the table for the next potential
        param_tab.clear();

        if(_genref){ // write ref aa ensemble cg-cg histograms

            file_name = (*potiter)->potentialName;
            file_name = file_name + aahistfile_extension;

            // print output file names on stdout
            cout << "Writing file: " << file_name << endl;
            
            // copy table from potential histogram
            aahist_tab = (*potiter)->aahist.data();

            // save table in the file
            aahist_tab.Save(file_name);
            // clear the table for the next potential
            aahist_tab.clear();

        }
        
    }
    
}

void CsgREupdate::EvalConfiguration(Topology *conf, Topology *conf_atom){

    /* as per Relative Entropy Ref.  J. Chem. Phys. 134, 094112, 2011
     * for each CG we need to compute dU/dlamda and d2U/dlamda_i dlamda_j
     * here, we add running sum of the second term of eq. 51 and store it
     * in _DS, ensemble avergaging and first term addition is done in EndEvaluate
     * for eq. 52, here the running sum of second and third terms is computed
     * and added to _HS, addition of 1st and 4th term and ensemble average
     * is performed in EndEvalute
     */
    // so lets do it

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

     if(!_genref) { // if not generating ref aa hist
        // update _DS and _HS
        for( int row = 0; row < _nlamda; row++ ){

            _DS(row) += (-1.0 * _beta * _dUFrame(row));

            for( int col = 0; col < _nlamda; col++){

                _HS(row,col) += (_beta * _beta * _dUFrame(row) * _dUFrame(col));

            }
        }

     } else { // map the first frame to cg topology
              // this topology is then used as the starting config for cg runs
         if(_nframes == 1)
             _writer->Write(conf);
         
     }
     
     _nframes++;
}

CsgREupdate::PotentialInfo::PotentialInfo(Topology *top, int index,
                                          bool bonded_, bool genref,
                                          int vec_pos_, Property* options) {
    //cout << "Potential info constructor" << endl;
    potentialIndex = index;
    _options = options;
    bonded   = bonded_;
    vec_pos  = vec_pos_;
    potentialName = options->get("name").value();
    type1 = options->get("type1").value();
    type2 = options->get("type2").value();
    potentialFunction = options->get("function").value();

    // count no. of beads
    BeadList beads1, beads2;
    beads1.Generate(*top, type1);
    beads2.Generate(*top, type2);
    // check beads exists
    if(beads1.size() == 0)
            throw std::runtime_error("Topology does not have beads of type \""
                        + type1 + "\"\n"
                    "This was specified in type1 of interaction \""
                    + potentialName + "\"");
    if(beads2.size() == 0)
            throw std::runtime_error("Topology does not have beads of type \""
                        + type2 + "\"\n"
                    "This was specified in type2 of interaction \""
                    + potentialName + "\"");

    nbeads1 = beads1.size();
    nbeads2 = beads2.size();

    // compute output table grid points
    double rmin = options->get("min").as<double>();
    double rmax = options->get("max").as<double>();
    double dr   = options->get("step").as<double>();
    int ngrid   = (int)( (rmax - rmin)/dr + 1.00000001 );

    pottblgrid.resize(ngrid,false);

    double r_init;
    int i;
    for (r_init = rmin, i = 0; i < ngrid - 1; r_init += dr) {

        pottblgrid(i++) = r_init;

    }
    pottblgrid(i) = rmax;
    // assign the user selected function form for this potential
    if( potentialFunction == "LJ126")
        ucg = new FunctionLJ126(rmin, rmax);
    else if (potentialFunction == "LJG")
        ucg = new FunctionLJG(rmin, rmax);
    else if (potentialFunction == "BSPL"){
        
        // get number of B-splines coefficients which are to be optimized
        int nlam = options->get("n_bspl_coeff").as<int>();
        // get number of B-splines coefficients near cut-off which must be
        // fixed to zero to ensure potential and forces smoothly go to zero
        int ncut = options->get("n_bspl_cut_coeff").as<int>();

        ucg = new FunctionBSPL(nlam, ncut, rmin, rmax);

    }
    else
        throw std::runtime_error("Function form \""
                        + potentialFunction + "\" selected for \""
                    + potentialName + "\" is not available yet.\n"
                    + "Please specify either \"LJ126, LJG, LJ2G, LJGC7, or BSPL\" "
                    + "in options file.");

    // initialize cg potential with old parameters
    
    string oldparam_file_extension = ".param.cur";
    string oldparam_file_name = potentialName + oldparam_file_extension;

    // if function form is BSPL and this is step 0 it must be fitted to given potential
    if( genref && potentialFunction == "BSPL" ){

        ucg->fitParam(oldparam_file_name);

    } else {

        ucg->setParam(oldparam_file_name);
        
    }
    
    /* read/load histogram if doing csg_reupdate
     * else program generates the histogram
     */
    string aahist_file_extension = ".aa.nbhist";
    string aahist_file_name = potentialName + aahist_file_extension;
 
    if(genref)
        // aahist uses the same rmin,rmax, and grid points as pot tables
        aahist.Initialize(rmin,rmax,ngrid);
    else
        aahist.Initialize(aahist_file_name);

    // let user know the properties of CG potential he/she selected
    cout << "We have " << potentialName << " CG potential" << endl;
    cout << "\t \t Between beads " << type1 << "-" << type2 << endl;
    cout << "\t \t With Function form " << potentialFunction << endl;
    cout << "\t \t And " << ucg->getParamSize() << " parameters" << endl;
    cout << "\t \t " << type1 << " has " << nbeads1 << " beads" << endl;
    if ( type1 != type2)
        cout << "\t \t " << type2 << " has " << nbeads2 << " beads" << endl;
    cout << "Potential range:" << endl;
    cout << "\t \t rmin    = " << rmin << " [nm]" << endl;
    cout << "\t \t rcutoff = " << rmax << " [nm]" << endl;
    cout << "\t \t step    = " << dr   << " [nm]" << endl;

}

// load use provided .xml option file
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
     
     /*
     for( int row = 0; row < _nlamda; row++){
         for(int col = 0; col < _nlamda; col++){
             cout << _HS(row,col)  << "\t" << _HS(col,row) << endl;
         }

     }
      */
     

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
    catch (NOT_SYM_POS_DEF){
        /* then can not use Newton-Raphson
         * go along steepest descent i.e. _dlamda = -_DS
         */
        cout << "**** _HS NOT positive definite ****" << endl;
        cout << "**** Using steepest descent! ****" << endl;
        
        _dlamda = minusDS;
    }
    // check if user opted to do convergence check
    if(_options.get("cg.inverse.convergence_check.type").as<string>() == "re"){

        CheckConvergence();

    }
    _lamda = _lamda + _relax * _dlamda ;

    // now update parameters of individual cg potentials
    PotentialContainer::iterator potiter;

    for (potiter = _potentials.begin();
            potiter != _potentials.end(); ++potiter){

        int pos_start = (*potiter)->vec_pos;
        int pos_max  = pos_start + (*potiter)->ucg->getParamSize();

        for( int row = pos_start; row < pos_max; row++){

            int lamda_i = row - pos_start;
            (*potiter)->ucg->setParam(lamda_i,_lamda(row));
            
        } // end row loop

    }// end potiter loop
    
}
//do nonbonded potential related update stuff for the current frame in evalconfig
// need to correct this routine for mixtures..seems ok for single fluid
void CsgREupdate::EvalNonbonded(Topology* conf, PotentialInfo* potinfo) {

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

    if (!_genref){ // if not generating ref aa histograms
        int pos_start = potinfo->vec_pos;
        int pos_max   = pos_start + potinfo->ucg->getParamSize();
        int lamda_i,lamda_j;
        double dU_i, d2U_ij;
        double U;

        // compute total energy
        U = 0.0;
        for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {
                U += potinfo->ucg->CalculateF((*pair_iter)->dist());
            }
        _UavgCG += U;
    
        // computing dU/dlamda and d2U/dlamda_i dlamda_j
        for( int row = pos_start; row < pos_max; row++){
        
            // ith parameter of this potential
            lamda_i = row - pos_start;

            // compute dU/dlamda
            // double sum by iterating over all pairs
            dU_i = 0.0;
            for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {

                dU_i += potinfo->ucg->CalculateDF(lamda_i,(*pair_iter)->dist());

            } // end loop pair_iter
        
            _dUFrame(row) = dU_i;
        
            for( int col = pos_start; col < pos_max; col++){
            
                lamda_j = col - pos_start;
      
                // compute d2U/dlamda_i dlamda_j and add to _HS
                // double sum by iterating over all pairs
                d2U_ij = 0.0;
                for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {

                    d2U_ij += potinfo->ucg->CalculateD2F(lamda_i,lamda_j,(*pair_iter)->dist());

                } // end loop pair_iter

                _HS(row,col) += (-1.0*_beta*d2U_ij);

            } // end loop col

        } // end loop row

    } else { // gather aa histogram data

        for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {
                potinfo->aahist.Process((*pair_iter)->dist(),1);
            }
    }

    
    delete nb;
}

// do non bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgNonbonded(PotentialInfo* potinfo) {

    /* here we need only number of type1 beads
     * since histogram is for typ1-type2
     * so nbeads = number of type1 beads
     * if both are same divide nbeads by 2 to avoid duplication
     */
    
    int pos_start = potinfo->vec_pos;
    int pos_max   = pos_start + potinfo->ucg->getParamSize();
    int lamda_i,lamda_j;
    double dU_i,d2U_ij;
    double U;
    Table hist;
    hist = potinfo->aahist.data();
    // compute avg AA energy
    U = 0.0;
    for(int bin = 0; bin < hist.size(); bin++) {

        double r_hist = hist.x(bin);
        double n_hist = hist.y(bin);
        U +=  n_hist *  potinfo->ucg->CalculateF(r_hist);

    }
    _UavgAA += U;
    
    // computing dU/dlamda and d2U/dlamda_i dlamda_j
    for( int row = pos_start; row < pos_max; row++){

        // ith parameter of this potential
        lamda_i = row - pos_start;

        // compute dU/dlamda and add to _DS
        dU_i = 0.0;
        for(int bin = 0; bin < hist.size(); bin++) {

            double r_hist = hist.x(bin);
            double n_hist = hist.y(bin);
            dU_i += n_hist * potinfo->ucg->CalculateDF(lamda_i,r_hist);
            
        } // end loop over hist
        _DS(row) += _beta * dU_i;
        
        for( int col = pos_start; col < pos_max; col++){

            lamda_j = col - pos_start;
            
            // compute d2U/dlamda_i dlamda_j and add to _HS
            d2U_ij = 0.0;
            for(int bin = 0; bin < hist.size(); bin++) {

                double r_hist = hist.x(bin);
                double n_hist = hist.y(bin);
                d2U_ij +=  n_hist *
                                potinfo->ucg->CalculateD2F(lamda_i,lamda_j,r_hist);

            } // end loop pair_iter

            _HS(row,col) += _beta * d2U_ij;
            
        } // end loop col

    } // end loop row
    
}

//do bonded potential related update stuff for the current frame in evalconfig
void CsgREupdate::EvalBonded(Topology* conf, PotentialInfo* potinfo){

    
}

// do bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgBonded(PotentialInfo* potinfo) {

}

// Compute refence AA ensemble CG-CG pair neibhor distances histogram
void CsgREupdate::AAavgHist(){

     // loop over all CG pairs to write generate aa ensemble histogram
     PotentialContainer::iterator potiter;

     for (potiter = _potentials.begin();
             potiter != _potentials.end(); ++potiter){

         PotentialInfo *potinfo = *potiter;
         potinfo->aahist.Normalize(_nframes);
     }

}

// Checks whether solution is converged using relative error 2-norm
void CsgREupdate::CheckConvergence(){

    /* if 2-norm of the relative error i.e. _dlamda/_lamda is less
     * than user specified tolerance value
     * program creates file "stop" in the main directory
     * which then is detected by inverse script and iterations stop
     */

    ub::vector<double> err;
    err.resize(_nlamda,false);
    err.clear();

    for( int row = 0; row < _nlamda; row++){

        if( _lamda(row) != 0.0 )
            err(row) = _dlamda(row)/_lamda(row);
        else
            err(row) = _dlamda(row);
    }

    double errnorm = ub::norm_inf(err);
    double tol = _options.get("cg.inverse.convergence_check.tol").as<double>();

    if( errnorm < tol ){
        
        ofstream out;
        string filename = "converged";
        out.open(filename.c_str());

        if(!out)
            throw runtime_error(string("error, cannot open file ") + filename);
        
        out << "Relative Entropy Converged!";

        out.close();
    }



}
