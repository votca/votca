/* 
 * File:   analysistool.h
 * Author: ruehle
 *
 * Created on August 2, 2007, 3:12 PM
 */

#ifndef _analasystool_H
#define	_analasystool_H

#include <map>
#include <string>
#include "cgengine.h"
#include "bondedstatistics.h"

using namespace std;

/**
    \brief base class for all analasys tools

    This is the base class for all analasys tool. 
    \todo do option functions!!!
*/
class AnalysisTool
{
public:
    AnalysisTool() {}
    virtual ~AnalysisTool() {}
    
    virtual void Register(map<string, AnalysisTool *> &lib) {}
    virtual void Command(BondedStatistics &bs, string cmd, vector<string> &args) {};
    virtual void Help(string cmd, vector<string> &args) {};
    
private:
//    map<string, string> _options;
};

#endif	/* _analasystool_H */

