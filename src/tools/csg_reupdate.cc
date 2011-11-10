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
//#include "doublesum_f.h"
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
    cout << "Initializing csg_reupdate" << endl;
     CsgApplication::Initialize();
    // add RE options
    AddProgramOptions("RE Specific options")
    ("options", boost::program_options::value<string>(), "  options file for coarse graining");
    
}

bool CsgREupdate::EvaluateOptions() {
    cout << "Evaluating options" << endl;
    CsgApplication::EvaluateOptions();
    CheckRequired("options", "need to specify options file");
    CheckRequired("trj", "no trajectory file specified");
    LoadOptions(OptionsMap()["options"].as<string>());

    return true;
}

void CsgREupdate::BeginEvaluate(Topology *top, Topology *top_atom){

    cout << "Beginning evaluate" << endl;

     _nlamda = 0;
    // initializing non-bonded interactions
     for (list<Property*>::iterator iter = _nonbonded.begin();
             iter != _nonbonded.end(); ++iter) {

         PotentialInfo *i = new PotentialInfo(_potentials.size(), false, _nlamda, *iter);
         // update parameter counter
         _nlamda += i->ucg.getParamSize();
         // add potential to container
         _potentials.push_back(i);
         
     }

     cout << "total parameters " << _nlamda << endl;

     _lamda.resize(_nlamda,false);
     // need to store initial guess of parameters in _lamda
    PotentialContainer::iterator potiter;

    for (potiter = _potentials.begin();
            potiter != _potentials.end(); ++potiter){

        int pos_start = (*potiter)->vec_pos;
        int pos_max  = pos_start + (*potiter)->ucg.getParamSize();

        for( int row = pos_start; row < pos_max; row++){

            int lamda_i = row - pos_start;
            _lamda(row) = (*potiter)->ucg.getParam(lamda_i);

        } // end row loop

    }// end potiter loop

    _DS.resize(_nlamda,false);
    _DS.clear();
    _HS.resize(_nlamda,_nlamda,false);
    _HS.clear();

    _nframes = 0.0; // no frames processed yet!
     // set Temperature
    _beta = _options.get("cg.inverse.kBT").as<double>();
     _relax = _options.get("cg.inverse.re.relax").as<double>();
     
     cout << "T in kBT " << _beta << endl;
}

void CsgREupdate::EndEvaluate(){
    cout << "Ending evaluate" << endl;

    // formulate _HS x = -_DS
    REFormulateLinEq();

    // solve _HS x = -_DS and update _lamda
    REUpdateLamda();
    
    WriteOutFiles();
    
    cout <<"Finished RE update!\n";
}

void CsgREupdate::WriteOutFiles() {

    cout << "Processed total of " << _nframes << " frames" << endl;
    cout<< "Writing updated CG potential table\n";

    string potfile_extension = ".pot.new";
    string paramfile_extension = ".param.new";
    string file_name;
    Table pot_tab;
    Table param_tab;
    // table with no error column
    pot_tab.SetHasYErr(false);
    param_tab.SetHasYErr(false);
    
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
            double out_u = (*potiter)->ucg.CalculateF(out_x);
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
        int nparam = (*potiter)->ucg.getParamSize();
        param_tab.resize(nparam, false);

        // print output file names on stdout
        cout << "Writing file: " << file_name << endl;

        // loop over output grid points
        for (int i = 0; i < nparam; i++) {

            int out_indx  = i;
            double out_param = (*potiter)->ucg.getParam(out_indx);
            // put point, result, flag at point out_x into the table
            param_tab.set(i, out_indx, out_param, 'i');

        }
        // save table in the file
        param_tab.Save(file_name);
        // clear the table for the next potential
        param_tab.clear();

    }
    
}

void CsgREupdate::EvalConfiguration(Topology *conf, Topology *conf_atom){

    
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

     _nframes++;
}

CsgREupdate::PotentialInfo::PotentialInfo(int index, 
                               bool bonded_, int vec_pos_, Property* options) {
    cout << "Potential info constructor" << endl;
    potentialIndex = index;
    _options = options;
    bonded   = bonded_;
    vec_pos  = vec_pos_;
    potentialName = options->get("name").value();
    type1 = options->get("type1").value();
    type2 = options->get("type2").value();

    cout << "type1 " << type1 << endl;
    cout << "type2 " << type2 << endl;
    
    string oldparam_file_extension = ".param.old";
    string oldparam_file_name = potentialName + oldparam_file_extension;

    // initialize cg potential with old parameters
    ucg.setParam(oldparam_file_name);

    // set min dist of potential
    ucg.setMinDist(options->get("min").as<double>());
    
    // set cut-off of potential
    ucg.setCutOffDist(options->get("cutoff").as<double>());

    string aahist_file_extension = ".aa.hist";
    string aahist_file_name = potentialName + aahist_file_extension;

    // read/load histogram
    aahist.Load(aahist_file_name);

    // compute output table grid points
    double rmin = options->get("pottbl.min").as<double>();
    double rmax = options->get("pottbl.max").as<double>();
    double dr   = options->get("pottbl.step").as<double>();
    int ngrid   = (int)( (rmax - rmin)/dr + 1.00000001 );

    pottblgrid.resize(ngrid,false);

    double r_init;
    int i;
    for (r_init = rmin, i = 0; i < ngrid - 1; r_init += dr) {

        pottblgrid(i++) = r_init;
        
    }
    pottblgrid(i) = rmax;

}

// load use provided .xml option file
void CsgREupdate::LoadOptions(const string &file) {
    cout << "Loading options" << endl;
    load_property_from_xml(_options, file);
    _nonbonded = _options.Select("cg.non-bonded");
}

// formulate _HS x = -_DS
void CsgREupdate::REFormulateLinEq() {

    /* compute CG ensemble avges of dU,d2U by dividing its
     * sum over all frames by total no. of frames
     */
    _DS /= (double)_nframes;
    _HS /= (double)_nframes;

    /* adding 4th term in eq. 52 of ref J. Chem. Phys. 134, 094112, 2011
     * to _HS
     */
    for( int row = 0; row < _nlamda; row++) {

        for( int col = 0; col < _nlamda; col++){

            _HS(row,col) += -1.0*_DS(row)*_DS(col);
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
    ub::vector<double> dlamda;
    dlamda.resize(_nlamda,false);
    ub::vector<double> residual(_nlamda);

    _DS *= -1.0;
    dlamda.clear();
    
    votca::tools::linalg_qrsolve(dlamda, _HS, _DS, &residual);

    _lamda += _relax * dlamda;

    // now update parameters of individual cg potentials
    PotentialContainer::iterator potiter;

    for (potiter = _potentials.begin();
            potiter != _potentials.end(); ++potiter){

        int pos_start = (*potiter)->vec_pos;
        int pos_max  = pos_start + (*potiter)->ucg.getParamSize();

        for( int row = pos_start; row < pos_max; row++){

            int lamda_i = row - pos_start;
            (*potiter)->ucg.setParam(lamda_i,_lamda(row));

        } // end row loop

    }// end potiter loop
    
}
//do nonbonded potential related update stuff for the current frame in evalconfig
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

    nb->setCutoff(potinfo->ucg.getCutOff());

    if( potinfo->type1 == potinfo->type2 ) { // same beads
        nb->Generate(beads1,true);
    } else { // different beads
        nb->Generate(beads1,beads2,true);
    }
    
    int pos_start = potinfo->vec_pos;
    int pos_max   = pos_start + potinfo->ucg.getParamSize();
    NBList::iterator pair_iter;
    int lamda_i,lamda_j;
    double dU_i,dU_j,d2U_ij;
  
    /* save some memory!
     * for eq 51 and 52 in Ref:: J. Chem. Phys. 134, 094112, 2011
     * now just compute running sum of the second term in eq. 51 and store in _DS
     * also compute running sum of (2nd + 3rd) term in eq. 52 and store in _HS
     * after all frame computations are finished rest of the terms in eq. 51,52
     * will be computed.
     */

    // computing dU/dlamda and d2U/dlamda_i dlamda_j
    for( int row = pos_start; row < pos_max; row++){

        dU_i = 0.0; dU_j = 0.0; d2U_ij = 0.0;
        
        // ith parameter of this potential
        lamda_i = row - pos_start;

        // compute dU/dlamda and add to _DS
        // double sum by iterating over all pairs
        for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {

            // assuming that nb->dist() takes care of pbc
            dU_i += potinfo->ucg.CalculateDF(lamda_i,(*pair_iter)->dist());
            
        } // end loop pair_iter

        _DS(row) += -1.0*_beta*dU_i;
        
        for( int col = pos_start; col < pos_max; col++){

            lamda_j = col - pos_start;
      
            // compute d2U/dlamda_i dlamda_j and add to _HS
            // double sum by iterating over all pairs
            for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {

                d2U_ij += potinfo->ucg.CalculateD2F(lamda_i,lamda_j,(*pair_iter)->dist());
                dU_j   += potinfo->ucg.CalculateDF(lamda_j,(*pair_iter)->dist());

            } // end loop pair_iter

            _HS(row,col) += -1.0*_beta*d2U_ij + _beta*_beta*dU_i*dU_j;

        } // end loop col

    } // end loop row

    delete nb;
}

// do non bonded potential AA ensemble avg energy computations
void CsgREupdate::AAavgNonbonded(PotentialInfo* potinfo) {

    int pos_start = potinfo->vec_pos;
    int pos_max   = pos_start + potinfo->ucg.getParamSize();
    int lamda_i,lamda_j;
    double dU_i,d2U_ij;
    
    // computing dU/dlamda and d2U/dlamda_i dlamda_j
    for( int row = pos_start; row < pos_max; row++){

        dU_i = 0.0; d2U_ij = 0.0;
        // ith parameter of this potential
        lamda_i = row - pos_start;

        // compute dU/dlamda and add to _DS
        for(int bin = 0; bin < potinfo->aahist.size(); bin++) {

            double r_hist = potinfo->aahist.x(bin);
            double n_hist = potinfo->aahist.y(bin);
            dU_i += n_hist * potinfo->ucg.CalculateDF(lamda_i,r_hist);

        } // end loop over hist
        
        _DS(row) += _beta * dU_i;

        for( int col = pos_start; col < pos_max; col++){

            lamda_j = col - pos_start;

            // compute d2U/dlamda_i dlamda_j and add to _HS
            for(int bin = 0; bin < potinfo->aahist.size(); bin++) {

                double r_hist = potinfo->aahist.x(bin);
                double n_hist = potinfo->aahist.y(bin);
                d2U_ij +=  n_hist *
                                potinfo->ucg.CalculateD2F(lamda_i,lamda_j,r_hist);

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
