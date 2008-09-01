/// \addtogroup csg
///@{
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
#include <tools/objectfactory.h>
#include "cgengine.h"
#include "bondedstatistics.h"

using namespace std;

/**
    \brief base class for all analasys tools

    This is the base class for all analasys tool. To implement a new analasys tool, derive it from this class
    and call RegisterObject(AnalasysFactory, UserTool, key) where UserTool is the name of the class....

    \todo do option functions!!!
*/
class AnalysisTool
{
public:
    AnalysisTool() {}
    virtual ~AnalysisTool() {}
    
    virtual void RegisteredAt(ObjectFactory<string, AnalysisTool> &factory) {}    
    virtual void Command(BondedStatistics &bs, string cmd, vector<string> &args) {};
    
private:
//    map<string, string> _options;
};

// important - singleton pattern, make sure factory is created before accessed
inline ObjectFactory<string, AnalysisTool> &AnalysisFactory()
{
    static ObjectFactory<string, AnalysisTool> _AnalysisFactory;
    return _AnalysisFactory;
}

#endif	/* _analasystool_H */

/// @}
