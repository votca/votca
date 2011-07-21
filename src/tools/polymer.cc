#include "polymer.h"

void Polymer::Initialize(){
    QMApplication::Initialize();
    /// define standard program options
    AddProgramOptions("Polymer options")
    ("cutnrg", boost::program_options::value<double>(), "wavefunctions will be considered only if their eigenvalues are greater than the smallest one by this amount");
 
}

bool Polymer::EvaluateOptions(){
    QMApplication::EvaluateOptions();
    CheckRequired("cutnrg", "no energy cutoff given");
    _cutnrg = OptionsMap()["cutnrg"].as<double>();
}

bool Polymer::EvaluateFrame(){
    UpdatePolTop();
    string outn ("frame");
    outn = outn + lexical_cast<string> (_qmtop.getStep());
    Save(outn);
    Clean();
}

void Polymer::EndEvaluate(){
    
}

void Polymer::UpdatePolTop(){
    // 1) compute "wavefunction" for each molecule 
    MoleculeContainer::iterator itMol=_qmtop.Molecules().begin();
    
    for ( ; itMol!= _qmtop.Molecules().end(); ++itMol){
        CalcWaveFunction(*itMol);
    }

    
    // 2) create new neighbour list
    QMNBList &nblist = _qmtop.nblist();
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter){
        CrgUnit * one = (*iter)->Crg1();
        CrgUnit * two = (*iter)->Crg2();
        //insert check for size of Js
        double J = ((*iter)->Js())[0];
        UpdateJs (one, two,J); // hard code that
                               // only one element is 
                               // considered (otherwise
                               // the whole method is doubtful)
    }

    // 3) now that the new neighbour list is there, compute the dr.
    map<PairWF, double>::iterator itpwf=_polJs.begin();
    for(;itpwf != _polJs.end(); ++itpwf){
        UpdatedR(itpwf->first);
    }
}


void Polymer::UpdateJs(CrgUnit * one, CrgUnit *two, const double & J){

    int mol1 = one->getMolId();
    int mol2 = two->getMolId();
    if (mol1 == mol2) {
        return;
    }
    vector <WaveFunction *>::iterator it1=_mstates[mol1]->begin();
    vector <WaveFunction *>::iterator it2=_mstates[mol2]->begin();
    map<int, int>::iterator itm1;
    map<int, int>::iterator itm2;
    itm1=_mcrg2bs.find(one->getId());
    itm2=_mcrg2bs.find(two->getId());
    if (itm1 == _mcrg2bs.end() || itm2 == _mcrg2bs.end() ){
        throw runtime_error ("could not find the crgunit in the index in UpdateJs");
    }
    for ( ; it1 != _mstates[mol1]->end(); ++it1){
        for ( ; it2 != _mstates[mol2]->end(); ++it2){

            double dJ = gsl_vector_get( (*it1)->_wf,  itm1->second) *
                        gsl_vector_get( (*it2)->_wf,  itm2->second) * J;
            map < PairWF , double >::iterator
              itJ = _polJs.find( make_pair(*it1, *it2));
            if (itJ == _polJs.end()){
                _polJs.insert(make_pair (make_pair(*it1, *it2),dJ));
            }
            else {
                itJ->second += dJ;
            }

        }
    }
}

void Polymer::CalcWaveFunction(Molecule * mol){
    int nbeads = mol->BeadCount();
    static int id=0;
    gsl_matrix * pH = gsl_matrix_calloc( nbeads, nbeads);
    double Jmax = 0.5;
    vector <WaveFunction *> * melement = new vector < WaveFunction *> ;
    _mstates.push_back( melement);
    
    //1 set diagonal elements
    for (int i=0; i< nbeads;i++){
        QMCrgUnit * crg = dynamic_cast<QMBead*>(mol->getBead(i))
            ->GetCrgUnit();
		double nrg = crg->getTotalEnergy();
        gsl_matrix_set (pH,i,i, nrg);
        // set a map from CrgUnits -> the integer associated with them
        map<int, int>::iterator itm = _mcrg2bs.find(crg->getId());
        if (itm == _mcrg2bs.end()){
            _mcrg2bs.insert(make_pair(crg->getId(), i));
        }
        else {
            throw runtime_error ("each crgunit should only appeare once in CalcWaveFunction");
        }
        QMBead* tbead = dynamic_cast<QMBead*>(mol->getBead(i));
        if (tbead->getiPos() != 0 ){
            throw runtime_error("each qmbead should have ipos correpsonding to its pos in the molecule");
        }
    }

    //2 set offdiagonal elements
    // in principle it would be nice to do this via the interactions
    // but this is not copied/written in statesaver
    
    //InteractionContainer::iterator itC=mol->Interactions().begin();
    //for ( ;itC != mol->Interactions.end(); itC++){
    //}
    for (int i=0;i<nbeads-1;i++){
        int j = i+1;
        vec one;
        vec two;
        one = mol->getBead(i)->getU();
        two = mol->getBead(j)->getU();
        // this will change using the reverse polish notation object.
        // dribble dribble
        double J = Jmax * (one * two);
        gsl_matrix_set(pH, i,j, J);
        gsl_matrix_set(pH, j,i, J);
    }

    gsl_eigen_symmv_workspace * gslwork = gsl_eigen_symmv_alloc(nbeads);
    gsl_vector *eval = gsl_vector_alloc (nbeads);
    gsl_matrix *evec = gsl_matrix_alloc (nbeads, nbeads);
    gsl_eigen_symmv(pH,eval, evec,gslwork);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    
    //for each evector make a new wavefunction
    for (int i=0;i<nbeads;i++){
        if ( gsl_vector_get(eval,i) < gsl_vector_get(eval,0) + _cutnrg ){
            gsl_vector * pv = gsl_vector_alloc (nbeads);
            gsl_vector_memcpy(pv, &(gsl_matrix_column(evec, i).vector) );
            
            WaveFunction * wf = new WaveFunction;
            wf->_molid  = mol->getId();
            wf->_wfid = i;
            wf->_nrg = gsl_vector_get(eval,i);
            wf->_id = id;
            wf->_wf = pv;
            _states.push_back(wf);
            melement->push_back(wf);
            //determine the position
            wf->_pos=vec(0.,0.,0.);
            for (int j=0;j<nbeads;j++){
                wf->_pos += (gsl_vector_get(pv,j) * gsl_vector_get(pv,j)) *
                        (mol->getBead(j)->getPos());
            }
            id++;
        }
    }

    gsl_matrix_free(pH);
    gsl_matrix_free(evec);
    gsl_eigen_symmv_free(gslwork);
    gsl_vector_free(eval);
    
}

void Polymer::Save(string & outn){
    string outnrg = outn + string ("nrg.txt");
    ofstream  out;
    out.open(outnrg.c_str());
    vector <WaveFunction *>::iterator its = _states.begin();
    for (; its!=_states.end();its++ ){
        out << "mol id: " << (*its)->_molid
                << " wf id: " << (*its)->_wfid
                << " wf abs id: " << (*its)->_id
                << " energy " << (*its)->_nrg << '\n';
    }
    out << endl;
    out.close();

    string outneigh = outn + string ("neigh.txt");
    out.open(outneigh.c_str());
    map < PairWF, double> ::iterator itJ = _polJs.begin();
    map < PairWF, vec> ::iterator itR = _poldR.begin();
    for (; itJ!= _polJs.end(); ++itJ, ++itR){
        out << " wf ids: " << ((itJ->first).first)->_id <<
                " " << ((itJ->first).second)->_id <<
                " J " << (itJ->second) <<
                " dR " << (itR->second) <<'\n';
    }
    out.close();

}

void Polymer::Clean(){
    vector <WaveFunction *>::iterator its = _states.begin();
    for (; its!=_states.end();its++ ){
        gsl_vector_free((*its)->_wf);
        delete *its;
    }
    _states.clear();
    vector < vector < WaveFunction * > *>::iterator itms = _mstates.begin();
    for ( ; itms != _mstates.end(); ++itms){
        (*itms)->clear();
    }
    _mstates.clear();
    _polJs.clear();
    _poldR.clear();
    _mcrg2bs.clear();
}

void Polymer::UpdatedR(const PairWF & pwf){
    
    vec r1 = (pwf.first)->_pos;
    vec r2 =(pwf.second)->_pos;
    vec d = _qmtop.BCShortestConnection(r1, r2);
    _poldR.insert(make_pair(pwf, d));
    
}

int main(int argc, char** argv) {

    Polymer pol_calc;
    return pol_calc.Exec(argc, argv);
}
