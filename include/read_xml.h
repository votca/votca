/* 
 * File:   read_xml.h
 * Author: vehoff
 *
 * Created on May 26, 2010, 11:55 AM
 */

#ifndef _READ_XML_H
#define	_READ_XML_H

#include <votca/tools/parsexml.h>
#include "qmcalculator.h"
#include "qmpair.h"

class ReadXML : public QMCalculator, public ParseXML{
public:
    ReadXML() {};
    ~ReadXML() {};

    bool EvaluateFrame(int nr, int nframes, QMTopology *top);

private:
    /// input file stream
    ifstream _in_int;

    void ParseRoot(const string &el, map<string, string> &attr);
    void ParsePair(const string &el, map<string, string> &attr, QMNBList &nblist);
};

#endif	/* _READ_XML_H */

