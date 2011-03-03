/* 
 * File:   qmcrgunit.h
 * Author: ruehle
 *
 * Created on February 7, 2011, 2:33 PM
 */

#ifndef __VOTCA_MD2QM_QMCRGUNIT_H
#define	__VOTCA_MD2QM_QMCRGUNIT_H

#include <votca/moo/crgunit.h>

using namespace std;

class QMCrgUnit : public CrgUnit
{
public:
    QMCrgUnit() : _occupation_probability(0) {}
    QMCrgUnit(vector <vec> positions, vector <vec> norms, vector <vec> planes,
            const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId);

    QMCrgUnit(const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId);

    // don't forget to overload these!
    //void copyCrgUnit(CrgUnit & acrg);
    //void copyCrgUnit(CrgUnit & acrg, const int & id);

    double getOccupationProbability() { return _occupation_probability; }
    void  setOccupationProbability(double prob) { _occupation_probability = prob; }

protected:
    double _occupation_probability;

    
};

inline QMCrgUnit::QMCrgUnit(vector <vec> positions, vector <vec> norms, vector <vec> planes,
            const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId) : CrgUnit(positions, norms, planes, id, type, molId), _occupation_probability(0)
{
    
}

inline QMCrgUnit::QMCrgUnit(const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId) : CrgUnit(id, type, molId), _occupation_probability(0)
{
    
}


#endif	/* QMCRGUNIT_H */

