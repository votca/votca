#include "easyjobserver.h"
#include <votca/csg/nblist.h>
#include <qmnblist.h>
#include <votca/tools/histogramnew.h>
#include <kmc/vertex.h>
#include <kmc/hoppers.h>
#include <kmc/kmc.h>

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
    ofstream out_cont; // for velocity averaging
    ofstream out_diff; // for diffusion
    out_cont.open("kmc_cont.res");
    out_diff.open("kmc_diff.res");
    if(out_cont!=0){
        out_cont << "# number of run, time [sec], average velocity [m/s], mobility [cm^2/sec], electric field [V/m], end to end vector [nm]" << endl;
    }
    if(out_diff!=0){
        out_diff << "# time [sec], displacement in x, y and z direction [m]" << endl;
    }
    out_cont.close();
    out_diff.close();
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
            (*iter)->setJs(Js);
         }
    }
    /// calculate & check the rates
    CalcRates(nblist);
    MakeRatesSIUnits(nblist);
    print_nbs_to_file(nblist);
    /// create KMC graph
    graph kmc_grid;
    make_kmc_graph(&kmc_grid,nblist);
    /// run continuous KMC
    _out_cont.open("kmc_cont.res");
    _out_diff.open("kmc_diff.res");
    if(_out_cont!=0){
        _out_cont << "# number of run, time [sec], average velocity [m/s], mobility [cm^2/sec], electric field [V/m], end to end vector [nm]" << endl;
    }
    if(_out_diff!=0){
        _out_diff << "# time [sec], displacement in x, y and z direction [m]" << endl;
    }
    _out_cont.close();
    _out_diff.close();
    Random::init( 14, 122, 472, 1912 );
    kmc_grid.setGeneratorsOnly();
    hoppers charges(&kmc_grid);
    KMCAlg kmc_alg(1E-8, 1E-15, 1.2, 1, &charges);
    cout << "Starting continuos KMC..." << endl;
    kmc_alg.kmcPBC(10, _out_cont, _out_diff);
    cout << "Finished continuous KMC." << endl;

    /// Testing the state saver class
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
        double rate_12, rate_21 = 0.0;
        double Jeff2 = (*iter)->calcJeff2();
        //cout << "Jeff = " << Jeff << endl;
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
        if(MatchNNnames(crg1, crg2))
        {
            double prefactor = 1.0;
            /// reorganization energy in eV as given in list_charges.xml
            double reorg = 0.5 * (crg1->GetType()->GetReorg()+crg2->GetType()->GetReorg());
            /// free energy difference due to electric field, i.e. E*r_ij
            double dG_field = -_E * ((*iter)->r()) * RA * Ang;
            /// free energy difference due to different energy levels of molecules
            double dG_en = crg2->GetNRG() - crg1->GetNRG();
            /// electrostatics are taken into account in qmtopology and are contained in NRG
            /// total free energy difference
            double dG = dG_field + dG_en;
            /// Marcus rate from first to second
            rate_12 = prefactor * sqrt(PI/(reorg * _kT)) * Jeff2 *
                exp (-(dG + reorg)*(dG + reorg)/(4*_kT*reorg));
            /// Marcus rate from second to first (dG_field -> -dG_field)
            dG = -dG_field + dG_en;
            rate_21 = prefactor * sqrt(PI/(reorg * _kT)) * Jeff2 *
                exp (-(dG + reorg)*(dG + reorg)/(4*_kT*reorg));
        }
        (*iter)->setRate12(rate_12);
        (*iter)->setRate21(rate_21);
    }
}

void EasyJObserver::MakeRatesSIUnits(QMNBList &nblist){
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter)
    {
        (*iter)->rate12() *= 1/hbar;
        (*iter)->rate21() *= 1/hbar;
    }
}

void EasyJObserver::print_nbs_to_file(QMNBList &nblist){
    ofstream out_nbl;
    out_nbl.open("nbl_votca.res");
    if(out_nbl!=0){
        out_nbl << "Neighbours, J(0), J_eff, rate, r_ij, abs(r_ij) [Bohr]" << endl;
        QMNBList::iterator iter;
        for ( iter  = nblist.begin(); iter != nblist.end() ; ++iter){
            out_nbl << "(" << (*iter)->first->GetId() << "," << (*iter)->second->GetId() << "): ";
            out_nbl << (*iter)->Js()[0] << " " << sqrt((*iter)->calcJeff2()) << " " << (*iter)->rate12() << " ";
            out_nbl << (*iter)->r().getX() << " " << (*iter)->r().getY() << " " << (*iter)->r().getZ() << " ";
            out_nbl << " " << abs((*iter)->r()) << endl;
        }
    }
    out_nbl.close();
}

void EasyJObserver::make_kmc_graph(graph *a, QMNBList &nblist){
    cout << "[make_kmc_graph]: Building KMC Graph...";
    /// assign constants
    a->SetField(_E);
    /// set vertices equal to centers of mass
    BeadContainer::iterator it;
    for(it=_qmtop->Beads().begin(); it!=_qmtop->Beads().end(); ++it){
        a->AddVertex((*it)->getPos(),_E); /// TO DO: remove necessity for E-field at this point
    }
    /// set edges, two edges 1->2 and 2->1 are created per neighboring pair
    for(QMNBList::iterator iter = nblist.begin(); iter!=nblist.end();++iter){
        a->AddEdge((*iter)->first->GetId(), (*iter)->second->GetId(), (*iter)->rate12(),(*iter)->r());
        a->AddEdge((*iter)->second->GetId(), (*iter)->first->GetId(), (*iter)->rate21(),-(*iter)->r());
    }
    cout << " Done." << endl;
}