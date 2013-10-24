#include <votca/ctp/polarseg.h>


namespace votca { namespace ctp {

    
    
PolarSeg::PolarSeg(int id, vector<APolarSite*> &psites) : _id(id) {    
    for (int i = 0; i < psites.size(); ++i) {
        push_back(psites[i]);
    }
    this->CalcPos();
}


PolarSeg::PolarSeg(PolarSeg *templ) {
    
    for (int i = 0; i < templ->size(); ++i) {
        APolarSite *newSite = new APolarSite((*templ)[i]);
        push_back(newSite);
    }
    this->_id = templ->_id;
    this->_pos = templ->_pos;    
}
    
    
PolarSeg::~PolarSeg() {
   vector<APolarSite*> ::iterator pit;
   for (pit = begin(); pit < end(); ++pit) {         
       delete *pit;
   }
   clear();
   vector<PolarNb*>::iterator nit;
   for (nit = _nbs.begin(); nit < _nbs.end(); ++nit) 
       delete *nit;
   _nbs.clear();
}


void PolarSeg::CalcPos() {    
    _pos = vec(0,0,0);    
    for (int i = 0; i < this->size(); ++i) {        
        _pos += (*this)[i]->getPos();        
    }
    if (this->size() > 0)
        _pos /= double(this->size());
}


double PolarSeg::CalcTotQ() {
    double Q = 0.0;
    for (int i = 0; i < this->size(); ++i) {
        Q += (*this)[i]->getQ00();
    }
    return Q;
}


void PolarSeg::Translate(const vec &shift) {    
    for (int i = 0; i < size(); ++i) {
        (*this)[i]->Translate(shift);
    }
    _pos += shift;
}


void PolarSeg::CalcIsCharged() {
    _is_charged = false;
    for (int i = 0; i < size(); ++i) {
        if ((*this)[i]->IsCharged()) _is_charged = true;
    }
    return;
}


void PolarSeg::CalcIsPolarizable() {
    _is_polarizable = false;
    for (int i = 0; i < size(); ++i) {
        if ((*this)[i]->IsPolarizable()) _is_polarizable = true;
    }
    return;
}


void PolarSeg::ClearPolarNbs() {
    vector<PolarNb*>::iterator nit;
    for (nit = _nbs.begin(); nit < _nbs.end(); ++nit) 
        delete *nit;
    _nbs.clear();
    return;
}


void PolarSeg::PrintPolarNbPDB(string outfile) {    
    FILE *out;
    out = fopen(outfile.c_str(),"w");
    PolarSeg::iterator pit;
    vector<PolarNb*>::iterator nit;
    for (pit = begin(); pit < end(); ++pit) {
        (*pit)->WritePdbLine(out, "CEN");
    }
    for (nit = _nbs.begin(); nit < _nbs.end(); ++nit) {
        PolarSeg *nb = (*nit)->getNb();
        nb->Translate((*nit)->getS());
        for (pit = nb->begin(); pit < nb->end(); ++pit) {
            (*pit)->WritePdbLine(out, "PNB");
        }
        nb->Translate(-1*(*nit)->getS());
    }
    fclose(out);
    return;
}


//template void PolarSeg::serialize<boost::archive::text_oarchive>
//    (boost::archive::text_oarchive &, const unsigned int);
//template void PolarSeg::serialize<boost::archive::text_iarchive>
//    (boost::archive::text_iarchive &, const unsigned int);
//
//template void PolarSeg::serialize<boost::archive::binary_oarchive>
//    (boost::archive::binary_oarchive &, const unsigned int);
//template void PolarSeg::serialize<boost::archive::binary_iarchive>
//    (boost::archive::binary_iarchive &, const unsigned int);


//template void PolarSeg::serialize<boost::archive::text_oarchive>
//    (boost::archive::text_oarchive &, const unsigned int);
//template void PolarSeg::serialize<boost::archive::text_iarchive>
//    (boost::archive::text_iarchive &, const unsigned int);
//
//template void PolarSeg::serialize<boost::archive::binary_oarchive>
//    (boost::archive::binary_oarchive &, const unsigned int);
//template void PolarSeg::serialize<boost::archive::binary_iarchive>
//    (boost::archive::binary_iarchive &, const unsigned int);

}}