#ifndef _WRITE_XML_H
#define	_WRITE_XML_H

#include "qmcalculator.h"
#include "qmpair.h"

/** \brief Writes out data to xml format

Callname: writexml

Writes out the following pair properties to xml format:
 - Distance between sites i and j in the pair in nm.
 - Distance vector coordinates x,y,z in nm.
 - Transfer integral in eV (or the effective value if more frontier orbitals are involved)
 - Intramolecular reorganization energy in eV.
 - Outer sphere reorganization energy in eV.

Writes out the following site properties to xml format:
 - Site energy in eV
 - Site occupation probability
 - Coordinates x,y,z of the center of mass in nm
*/
class WriteXML : public QMCalculator{
public:
    WriteXML() {};
    ~WriteXML() {};

    const char *Description() { return "Writes out data to xml format, TODO: remove this calculator once sqlite is running"; }

    void Initialize(QMTopology *top, Property *options);

    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);
    
private:
    string _outfile;
    bool _write_dist, _write_rij, _write_jeff, _write_en, _write_occ, _write_coords, _write_lambdaout;
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

    if(options->exists("options.writexml.occupation"))
        _write_occ = options->get("options.writexml.occupation").as<bool>();

    if(options->exists("options.writexml.positions"))
        _write_coords = options->get("options.writexml.positions").as<bool>();

    if(options->exists("options.writexml.lambdaout"))
        _write_lambdaout = options->get("options.writexml.lambdaout").as<bool>();


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
            << " first=\"" << (*iter)->first->getId() << "\""
            << " second=\"" << (*iter)->second->getId() << "\""
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
        if(_write_lambdaout)
            out << " lambdaout=\"" <<(*iter)->getLambdaOuter() << "\"";
        if(_write_jeff)
            out << " Jeff=\"" <<(*iter)->calcJeff2()<< "\"";
        out <<  "/>" << endl;
    }
        for (vector < QMCrgUnit *>::iterator lch_iter = lcharges.begin();lch_iter!=lcharges.end();++lch_iter) {
            out << "    <site "
                << " number=\"" << (*lch_iter)->getId() << "\"";
            if (_write_coords) {
                out << " pos=\"" << (*lch_iter)->GetCom() << "\"";
            }
            if (_write_en) {
                out << " energy=\"" << (*lch_iter)->getEnergy() << "\"";
            }
            if (_write_occ) {
                    out << " occupation=\"" << (*lch_iter)->getOccupationProbability() << "\"";
            }
            out << "/>" << endl;
    }
    out << "  </frame>" << endl;

    out.close();
}

#endif	/* _WRITE_XML_H */

