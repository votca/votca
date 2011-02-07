/* 
 * File:   write_xml.h
 * Author: vehoff
 *
 * Created on May 25, 2010, 3:01 PM
 */

#ifndef _WRITE_XML_H
#define	_WRITE_XML_H

#include "qmcalculator.h"
#include "qmpair.h"

class WriteXML : public QMCalculator{
public:
    WriteXML() {};
    ~WriteXML() {};

    void Initialize(QMTopology *top, Property *options);

    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);
    
private:
    string _outfile;
    bool _write_dist, _write_rij, _write_jeff, _write_en;
};

inline void WriteXML::Initialize(QMTopology *top, Property *options)
{
    _outfile = options->get("options.writexml.file").as<string>();
    _write_dist = false;

    if(options->exists("options.writexml.dist"))
        _write_dist = options->get("options.writexml.dist").as<bool>();

    if(options->exists("options.writexml.rij"))
        _write_rij = options->get("options.writexml.rij").as<bool>();

    if(options->exists("options.writexml.jeff"))
        _write_jeff = options->get("options.writexml.jeff").as<bool>();

    if(options->exists("options.writexml.site_energies"))
        _write_en = options->get("options.writexml.site_energies").as<bool>();


    ofstream out;
    out.open(_outfile.c_str(), ios::out);
    if(!out)
        throw std::runtime_error("error, cannot open " + _outfile);
    out << "<qmtop>" << endl;
    out.close();
}

inline void WriteXML::EndEvaluate(QMTopology *top)
{
    ofstream out;
    out.open(_outfile.c_str(), ios::out | ios::app);
    if(!out)
        throw std::runtime_error("error, cannot open " + _outfile);
    out << "</qmtop>" << endl;
    out.close();
}

inline bool WriteXML::EvaluateFrame(QMTopology *top){
    ofstream out;
    out.open(_outfile.c_str(), ios::out | ios::app);
    if(!out)
        throw std::runtime_error("error, cannot open " + _outfile);


    QMNBList &nblist = top->nblist();
    vector < QMCrgUnit *> lcharges = top->CrgUnits();
    
    out << "  <frame>"  << endl;
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter) {
        out << "    <pair "
            << " first=\"" << (*iter)->first->getId()+1 << "\""
            << " second=\"" << (*iter)->second->getId()+1 << "\""
            << " J=\"";
        for(int i=0; i< (*iter)->Js().size(); ++i)
            out << (*iter)->Js()[i] << " ";
        out << "\""
            << " rate12=\"" << (*iter)->rate12() << "\""
            << " rate21=\"" << (*iter)->rate21() << "\"";
        if(_write_dist)
            out << " dist=\"" <<(*iter)->dist()<< "\"";
        if(_write_rij)
            out << " rij=\"" <<(*iter)->r()<< "\"";
        if(_write_jeff)
            out << " Jeff=\"" <<(*iter)->calcJeff2()<< "\"";
        out <<  "/>" << endl;
    }
    if (_write_en) {
        for (vector < QMCrgUnit *>::iterator lch_iter = lcharges.begin();lch_iter!=lcharges.end();++lch_iter) {
            out << "    <site "
                << " number=\"" << (*lch_iter)->getId()+1 << "\""
                << " energy=\"" << (*lch_iter)->getEnergy() << "\"";
            out << "/>" << endl;
        }
    }
    out << "  </frame>" << endl;

    out.close();
}

#endif	/* _WRITE_XML_H */

