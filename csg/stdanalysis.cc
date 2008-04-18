// 
// File:   stdanalysis.cc
// Author: ruehle
//
// Created on August 2, 2007, 5:31 PM
//

#include <iostream>
#include <fstream>
#include <vector>
#include "analysistool.h"
#include <crosscorrelate.h>
#include <correlate.h>

class StdAnalysis
    : public AnalysisTool
{
    public:
        StdAnalysis() {};
        ~StdAnalysis() {};
        
        void RegisteredAt(ObjectFactory<string, AnalysisTool> &factory);
    
        void Command(CGEngine &cg, string cmd, vector<string> &args);                
        
        void WriteValues(CGEngine &cg, vector<string> &args);
        void WriteCorrelations(CGEngine &cg, vector<string> &args);
        void WriteAutocorrelation(CGEngine &cg, vector<string> &args);
    private:            
};

REGISTER_OBJECT(AnalysisFactory, StdAnalysis, "list");

void StdAnalysis::RegisteredAt(ObjectFactory<string, AnalysisTool> &factory)
{
    factory.Register("vals", this);
    factory.Register("cor", this);
    factory.Register("autocor", this);    
}

void StdAnalysis::Command(CGEngine &cg, string cmd, vector<string> &args)
{
    if(cmd == "vals") WriteValues(cg, args);
    if(cmd == "cor") WriteCorrelations(cg, args);
    if(cmd == "autocor") WriteAutocorrelation(cg, args);
    if(cmd == "list") {
        DataCollection<double>::selection *sel = cg.BondedValues().select("*");
        DataCollection<double>::selection::iterator i;
        cout << "Available bonded interactions:" << endl;
        for(i=sel->begin(); i!=sel->end(); ++i)
            cout << (*i)->getName() << " "; // << "[" << (*i).second->size() << "]" << " ";
        cout << endl;
        delete sel;
    }
}

void StdAnalysis::WriteValues(CGEngine &cg, vector<string> &args)
{
    ofstream out;
    
    DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = cg.BondedValues().select(args[i], sel);
    
    out.open(args[0].c_str());
    out << *sel << endl;
    out.close();
    cout << "written " << sel->size() << " data rows to " << args[0] << endl;
    delete sel;
}

void StdAnalysis::WriteAutocorrelation(CGEngine &cg, vector<string> &args)
{
    ofstream out;
    DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = cg.BondedValues().select(args[i], sel);
        
    CrossCorrelate c;
    c.AutoCorrelate(sel, false);
    out.open(args[0].c_str());
    out << c << endl;
    out.close();
    cout << "calculated autocorrelation for " << sel->size() << " data rows, written to " << args[0] << endl;
    delete sel;
}

void StdAnalysis::WriteCorrelations(CGEngine &cg, vector<string> &args)
{
    ofstream out;
        DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = cg.BondedValues().select(args[i], sel);
        
    Correlate c;
    c.CalcCorrelations(sel);
    out.open(args[0].c_str());
    out << c << endl;
    out.close();
    cout << "calculated correlations for " << sel->size() << " rows, written to " << args[0] << endl;
    delete sel;
}
