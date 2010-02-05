#include "integralsapp.h"

Integrals::Integrals()
{}

Integrals::~Integrals()
{}

void Integrals::HelpText()
{}

void Integrals::AddSpecificOptions(){
    namespace po = boost::program_options;

    /// define standard program options
    _op_desc_specific.add_options()
    ("nnnames", boost::program_options::value<string>()->default_value("*"), "  List of strings that the concatenation of the two molnames must match to be analyzed")
    ;
}

void Integrals::setNNnames(string nnnames){
    Tokenizer tok(nnnames, " ;");
    Tokenizer::iterator itok = tok.begin();
    for (; itok!= tok.end(); ++itok){
        _nnnames.push_back(*itok);
    }
}

void Integrals::Initialize(){
    /// set the nearest neighbor names to be matched
    setNNnames(_op_vm["nnnames"].as<string>());
}

bool Integrals::MatchNNnames(CrgUnit *crg1, CrgUnit* crg2){
    vector <string>::iterator its = _nnnames.begin();
    string namecrg = crg1->getType()->GetName()+string(":")+crg2->getType()->GetName();
    for ( ; its!= _nnnames.end(); ++its){
        if(wildcmp(its->c_str(), namecrg.c_str()) ){
            return true;
        }
    }
    return false;
}

bool Integrals::EvaluateFrame(){
    QMNBList &nblist = _qmtop.nblist();
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter){
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
        if(MatchNNnames(crg1, crg2)){
            vector <double> Js = _qmtop.GetJCalc().CalcJ(*crg1, *crg2);
            (*iter)->setJs(Js);
        }
    }
}

 void Integrals::EndEvaluate()
 {
     PrintNbs("nbl.res");
 }