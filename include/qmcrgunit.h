/* 
 * File:   qmcrgunit.h
 * Author: ruehle
 *
 * Created on February 7, 2011, 2:33 PM
 */

#ifndef __VOTCA_MD2QM_QMCRGUNIT_H
#define	__VOTCA_MD2QM_QMCRGUNIT_H

#include <votca/moo/crgunit.h>
#include <string>
#include <map>

using namespace std;

class QMCrgUnit : public CrgUnit
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
    
    // for Falk: clean out the old SetEnergy / getEnergy and then rename these functions
    void setEnergy_(string label, double value) { _energies[label] = value; }
    double getEnergy_(string label) { return _energies[label]; }
    map<string, double> &getEnergies() { return _energies; }
    
    double SumEnergies() { // TODO: add label + wildcard
        double e=0;
        for(map<string, double>::iterator iter = _energies.begin();
            iter!= _energies.end(); ++iter) e+= iter->second;
        return e;
    }

    // use pair<double, double> for lambdea. pair.first = charging, pair.second = discharging
    // setReorg(string label, double charging, double discharging) { _reorg[label] = makepair(charging, discharging); }
    // pair<double,double> &getReorg(string label) { return _reorg[label]; }
protected:
    double _occupation_probability;
    bool _in_database;

    map<string, double> _energies;
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


#endif	/* QMCRGUNIT_H */

