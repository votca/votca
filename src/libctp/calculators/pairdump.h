/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _PAIR_DUMP_H
#define	_PAIR_DUMP_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/paircalculator.h>
#include <sys/stat.h>
#include <votca/csg/trajectorywriter.h>

namespace votca { namespace ctp {

class PairDump : public PairCalculator
{
public:
    PairDump() {};
    ~PairDump() {};

    const char *Description() { return "Write pairs for DFT integral input"; }

    void EvaluateSite(QMTopology *top, QMCrgUnit *crg);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

   void Initialize(QMTopology *top, Property *options);

protected:
    TrajectoryWriter *_writer;
    string _format; // extension for writing files
    string _framedir;
    bool _subfolders;
    bool _mols;
    bool _pairs;
};

inline void PairDump::Initialize(QMTopology *top, Property *options)
{
    _format = "pdb";
    if(options->exists("options.pairdump.format"))
            _format = options->get("options.pairdump.format").as<string>();
    _subfolders=false;
    if(options->exists("options.pairdump.subfolders"))
            _subfolders= options->get("options.pairdump.subfolders").as<bool>();
    _mols=false;
    _pairs=true;
    if(options->exists("options.pairdump.molecules"))
            _mols = options->get("options.pairdump.molecules").as<bool>();
    if(options->exists("options.pairdump.pairs"))
            _pairs = options->get("options.pairdump.pairs").as<bool>();
        
}

inline bool PairDump::EvaluateFrame(QMTopology *top) {
    _writer = TrjWriterFactory().Create("." + _format);
    if(_writer == NULL)
        throw runtime_error(string("output format not supported: ") + _format);
    _framedir=string("frame")+boost::lexical_cast<string>(top->getStep()+1) +string("/") ;
    mkdir(_framedir.c_str(),0755);
    
    if(_pairs)
        PairCalculator::EvaluateFrame(top);

    if(!_mols) return true;
    vector<QMCrgUnit*>& crglist = top->CrgUnits();
    for(vector<QMCrgUnit*>::iterator iter = crglist.begin();iter!=crglist.end();++iter)
        EvaluateSite(top, *iter);
    return true;
}

inline void PairDump::EndEvaluate(QMTopology *top)
{
    delete _writer;
}

inline void PairDump::EvaluateSite(QMTopology *top, QMCrgUnit *crg)
{
    Topology atop;

    top->AddAtomisticBeads(crg,&atop);
    string subdir = _framedir;
    if(_subfolders) {
        subdir=subdir + "mol_" + boost::lexical_cast<string>(crg->getId()) + "/";
        mkdir(subdir.c_str(),0755);
    }

    string file = subdir
            + "mol_" + boost::lexical_cast<string>(crg->getId()) + "."  + _format;
    _writer->Open(file);
    _writer->Write(&atop);
    _writer->Close();
}

inline void PairDump::EvaluatePair(QMTopology *top, QMPair *pair){
    Topology atop;

    CrgUnit *crg1 = pair->Crg1PBCCopy();
    CrgUnit *crg2 = pair->Crg2PBCCopy();

    top->AddAtomisticBeads(crg1,&atop);
    top->AddAtomisticBeads(crg2,&atop);

    string subdir = _framedir;
    if(_subfolders) {
        subdir=subdir+
            + "pair_" + boost::lexical_cast<string>(crg1->getId()) + string("_")
            + boost::lexical_cast<string>(crg2->getId()) + "/";
        mkdir(subdir.c_str(),0755);
        subdir=subdir+"dim/";
        mkdir(subdir.c_str(),0755);
    }

    string file = subdir
            + "pair_" + boost::lexical_cast<string>(crg1->getId()) + string("_")
            + boost::lexical_cast<string>(crg2->getId()) + "."  + _format;
    _writer->Open(file);
    _writer->Write(&atop);
    _writer->Close();
}

}}

#endif	
