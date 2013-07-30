/* 
 * File:   qmpackage.h
 * Author: andrienko
 *
 * Created on July 30, 2013, 8:25 PM
 */

#ifndef _CTP_QMPACKAGE_H
#define	_CTP_QMPACKAGE_H

#include <votca/ctp/logger.h>

namespace votca { namespace ctp {
    
class QMPackage
{

public:

    virtual std::string getPackageName() = 0;
    virtual bool Run() = 0;
    virtual bool WriteInputFile() = 0;
    virtual bool ParseOutput() = 0;
    
    void setLog( Logger *pLog ) { _pLog = pLog; }
    Logger *getLog() { return _pLog; }
   
private:
    
    Logger *_pLog;
    
};

}}

#endif	/* QMPACKAGE_H */

