/* 
 * File:   integralsapp.h
 * Author: vehoff
 *
 * Created on February 3, 2010, 5:59 PM
 */

#ifndef _INTEGRALSAPP_H
#define	_INTEGRALSAPP_H

#include "qmapplication.h"
#include <sstream>

class Integrals : public QMApplication
{
public:
    Integrals();
    ~Integrals();

    void HelpText();
    void AddSpecificOptions();
    void Initialize();
    void EvaluateFrame();
    void EndEvaluate();
    
private:
    /// set _nnnames based on input
    void setNNnames(string  nnnames);
    /// check nearest neighbor names
    bool MatchNNnames(CrgUnit *crg1, CrgUnit* crg2);

};

#endif	/* _INTEGRALSAPP_H */

