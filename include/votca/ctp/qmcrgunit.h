/* 
 * File:   qmcrgunit.h
 * Author: ruehle
 *
 * Created on February 7, 2011, 2:33 PM
 */

#ifndef __VOTCA_MD2QM_QMCRGUNIT_H
#define	__VOTCA_MD2QM_QMCRGUNIT_H

#include"customfields.h"
#include <votca/moo/crgunit.h>
#include <string>
#include <map>

namespace votca { namespace ctp {

using namespace std;
using namespace votca::moo;

class QMCrgUnit : public CrgUnit, public CustomFields
{
public:
    QMCrgUnit() : _occupation_probability(0.), _in_database(false) {}
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

    void setInDatabase(bool indb) { _in_database = indb; }
    bool getInDatabase() { return _in_database; }       

    double getTotalEnergy(const char *match="energy_*", const char *nomatch="") {        
        //return getDouble("energy_intra")+getDouble("energy_coulomb");
        double e=0;
        for(map<string, double>::iterator iter = DoubleValues().begin();
            iter!= DoubleValues().end(); ++iter) {
            if(wildcmp(match, iter->first.c_str()))
                if(!wildcmp(nomatch, iter->first.c_str()))
                    e+= iter->second;
        }
        return e;
    }
    
protected:
    double _occupation_probability;
    bool _in_database;
};

inline QMCrgUnit::QMCrgUnit(vector <vec> positions, vector <vec> norms, vector <vec> planes,
            const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId) : CrgUnit(positions, norms, planes, id, type, molId), _occupation_probability(0.)
{
    
}

inline QMCrgUnit::QMCrgUnit(const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId) : CrgUnit(id, type, molId), _occupation_probability(0), _in_database(false)
{
    
}

}}

#endif	/* QMCRGUNIT_H */

