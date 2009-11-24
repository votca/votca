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
#include <votca/tools/crosscorrelate.h>
#include <votca/tools/correlate.h>
#include "bondedstatistics.h"
#include "stdanalysis.h"


void StdAnalysis::Register(map<string, AnalysisTool *> &lib)
{
    lib["list"] = this;
    lib["vals"] = this;
    lib["cor"] = this;
    lib["autocor"] = this;
}

void StdAnalysis::Command(BondedStatistics &bs, string cmd, vector<string> &args)
{
    if(cmd == "vals") WriteValues(bs, args);
    if(cmd == "cor") WriteCorrelations(bs, args);
    if(cmd == "autocor") WriteAutocorrelation(bs, args);
    if(cmd == "list") {
        DataCollection<double>::selection *sel = bs.BondedValues().select("*");
        DataCollection<double>::selection::iterator i;
        cout << "Available bonded interactions:" << endl;
        for(i=sel->begin(); i!=sel->end(); ++i)
            cout << (*i)->getName() << " "; // << "[" << (*i).second->size() << "]" << " ";
        cout << endl;
        delete sel;
    }
}

void StdAnalysis::WriteValues(BondedStatistics &bs, vector<string> &args)
{
    ofstream out;
    
    DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = bs.BondedValues().select(args[i], sel);
    
    out.open(args[0].c_str());
    out << *sel << endl;
    out.close();
    cout << "written " << sel->size() << " data rows to " << args[0] << endl;
    delete sel;
}

void StdAnalysis::WriteAutocorrelation(BondedStatistics &bs, vector<string> &args)
{
    ofstream out;
    DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = bs.BondedValues().select(args[i], sel);
        
    CrossCorrelate c;
    c.AutoCorrelate(sel, false);
    out.open(args[0].c_str());
    out << c << endl;
    out.close();
    cout << "calculated autocorrelation for " << sel->size() << " data rows, written to " << args[0] << endl;
    delete sel;
}

void StdAnalysis::WriteCorrelations(BondedStatistics &bs, vector<string> &args)
{
    ofstream out;
        DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = bs.BondedValues().select(args[i], sel);
        
    Correlate c;
    c.CalcCorrelations(sel);
    out.open(args[0].c_str());
    out << c << endl;
    out.close();
    cout << "calculated correlations for " << sel->size() << " rows, written to " << args[0] << endl;
    delete sel;
}
