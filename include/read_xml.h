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

class ReadXML : public QMCalculator
{
public:
    ReadXML() {};
    ~ReadXML() {};

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

private:
    void ParseRoot(const string &el, map<string, string> &attr);
    void ParseBody(const string &el, map<string, string> &attr);
    void ParseFrame(const string &el, map<string, string> &attr);
    void ParsePair(const string &el, map<string, string> &attr);

    ParseXML _parser;

    QMTopology *_top;
};

#endif	/* _READ_XML_H */

