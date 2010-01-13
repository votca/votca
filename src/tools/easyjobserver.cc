#include "easyjobserver.h"
#include <votca/csg/nblist.h>
#include <qmnblist.h>
#include <votca/tools/histogramnew.h>

EasyJObserver::EasyJObserver()
{}

EasyJObserver::~EasyJObserver()
{}

void EasyJObserver::Initialize(QMTopology &qmtop, Property &opts)
{
    _qmtop = &qmtop;
    _kT = opts.get("calc_rates.thermal_energy").as<double>();
    _E = opts.get("calc_rates.e_field").as<vec>();
    cout << "E = " << _E << endl;
}

void EasyJObserver::setNNnames(string  nnnames){
    Tokenizer tok(nnnames, " ;");
    Tokenizer::iterator itok = tok.begin();
    for (; itok!= tok.end(); ++itok){
        _nnnames.push_back(*itok);
    }
}
void EasyJObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);
}

void EasyJObserver::EndCG()
{
    HistogramNew histo;
    histo.Initialize(*(std::min_element( _Js.begin(), _Js.end() )),*(std::max_element( _Js.begin(), _Js.end() ) ) , int(_Js.size() / 10 ));
    histo.ProcessRange(_Js.begin(),_Js.end());
    cout << histo <<endl;
}

/// evaluate current conformation

void EasyJObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{
    _qmtop->Update(*top);
    QMNBList &nblist = _qmtop->nblist();

    BeadList list1;
    Topology *toptmp = dynamic_cast<Topology*>(_qmtop);
    list1.Generate(*toptmp, "*");

    nblist.setCutoff(_cutoff);
    nblist.Generate(list1);

    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
         if(MatchNNnames(crg1, crg2)){
            vector <double> Js = _qmtop->GetJCalc().GetJ(*crg1, *crg2);
            //cout << crg1->GetId() << " " << crg2->GetId() << " ";
            for(int i=0; i<Js.size(); ++i)
            {
                _Js.push_back(Js[i]);
            }
            (*iter)->setJ(Js[0]);
         }
    }
    CalcRates(nblist);
    MakeRatesSIUnits(nblist);
    cout<<"Falks test\n";
       StateSaver _saver(*_qmtop);
            string outfilename = "falks.dat";
            _saver.Open(outfilename,true);
  
          _saver.Write_QMBeads(_qmtop);
}

bool EasyJObserver::MatchNNnames(CrgUnit *crg1, CrgUnit* crg2){
    vector <string>::iterator its = _nnnames.begin();
    string namecrg = crg1->GetType()->GetName()+string(":")+crg2->GetType()->GetName();
    for ( ; its!= _nnnames.end(); ++its){
        if(wildcmp(its->c_str(), namecrg.c_str()) ){
            return true;
        }
    }
    return false;
}

void EasyJObserver::CalcRates(QMNBList &nblist){
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter)
    {
        double rate = 0.0;
        double Jeff = (*iter)->j();
        //cout << "Jeff = " << Jeff << endl;
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
        if(MatchNNnames(crg1, crg2))
        {
            double prefactor = 1.0;
            /// reorganization energy in eV as given in list_charges.xml
            double reorg = 0.5 * (crg1->GetType()->GetReorg()+crg2->GetType()->GetReorg());
            /// free energy difference due to electric field, i.e. E*r_ij
            double dG_field = _E * ((*iter)->r()) * RA * Ang;
            /// free energy difference due to different energy levels of molecules
            double dG_lev = crg2->GetNRG() - crg1->GetNRG();
            /// free energy difference due to electrostatics
            double dG_estatic = 0.0;
            /// total free energy difference
            double dG = dG_field + dG_lev + dG_estatic;
            /// Marcus rate
            rate = prefactor * sqrt(PI/(reorg * _kT) ) * Jeff*Jeff *
                exp (-(dG + reorg)*(dG + reorg)/(4*_kT*reorg));
            //cout << "Rate: " << rate << endl;
            //cout << "dG_field = " << dG_field << endl;
        }
        (*iter)->setRate(rate);
    }
}

void EasyJObserver::MakeRatesSIUnits(QMNBList &nblist){
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter)
    {
        (*iter)->rate() *= 1/hbar;
    }
}

