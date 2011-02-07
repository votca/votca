/* 
 * File:   calc_integrals.h
 * Author: vehoff
 *
 * Created on April 8, 2010, 10:33 AM
 */

#ifndef _PAIR_DUMP_H
#define	_PAIR_DUMP_H

#include "qmpair.h"
#include "paircalculator.h"
#include <sys/stat.h>
#include <votca/csg/trajectorywriter.h>

class PairDump : public PairCalculator
{
public:
    PairDump() {};
    ~PairDump() {};

    void EvaluatePair(QMTopology *top, QMPair *pair);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

   void Initialize(QMTopology *top, Property *options);

protected:
    TrajectoryWriter *_writer;
    string _format; // extension for writing files
    string _framedir;
};

inline void PairDump::Initialize(QMTopology *top, Property *options)
{
    _format = "pdb";
    if(options->exists("options..pairdump.format"))
            _format = options->get("options..pairdump.ext").as<string>();
}

bool PairDump::EvaluateFrame(QMTopology *top) {
    _writer = TrjWriterFactory().Create("." + _format);
    if(_writer == NULL)
        throw runtime_error(string("output format not supported: ") + _format);
    _framedir=string("frame")+lexical_cast<string>(top->getStep()) +string("/") ;
    mkdir(_framedir.c_str(),0755);
    PairCalculator::EvaluateFrame(top);
}

void PairDump::EndEvaluate(QMTopology *top)
{
    delete _writer;
}

inline void PairDump::EvaluatePair(QMTopology *top, QMPair *pair){
    Topology atop;

    CrgUnit *crg1 = pair->Crg1();
    CrgUnit *crg2 = pair->Crg2();

    top->AddAtomisticBeads(crg1,&atop);
    top->AddAtomisticBeads(crg2,&atop);

    string file = _framedir
            + boost::lexical_cast<string>(crg1->getId()+1) + string("_")
            + boost::lexical_cast<string>(crg2->getId()+1) + "."  + _format;
    _writer->Open(file);
    _writer->Write(&atop);
    _writer->Close();
}

#endif	
