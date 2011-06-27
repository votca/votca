#ifndef _READ_XML_H
#define	_READ_XML_H

#include <votca/tools/parsexml.h>
#include "qmcalculator.h"
#include "qmpair.h"

/** \brief Reads in data in xml format

Callname: readxml

Reads in the following pair properties to xml format:
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


class ReadXML : public QMCalculator
{
public:
    ReadXML() {};
    ~ReadXML() {};

    const char *Description() { return "Reads in data in xml format, TODO: remove this calculator once sqlite is running"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

private:
    void ParseRoot(const string &el, map<string, string> &attr);
    void ParseBody(const string &el, map<string, string> &attr);
    void ParseFrame(const string &el, map<string, string> &attr);
    void ParsePair(const string &el, map<string, string> &attr);

    ParseXML _parser;
    string _filename;
    
    QMTopology *_top;
};

#endif	/* _READ_XML_H */

