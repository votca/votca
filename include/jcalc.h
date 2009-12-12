/* 
 * File:   jcalc.h
 * Author: james
 *
 * Created on 07 December 2009, 17:19
 */

#ifndef _JCALC_H
#define	_JCALC_H

#include <boost/lexical_cast.hpp>
#include <libxml/parser.h>
#include <map>
#include <stdexcept>

#include "crgunittype.h"
#include "votca/tools/tokenizer.h"
#include "crgunit.h"
#include "fock.h"

using namespace boost;

class JCalc{

public:

    JCalc(){

    }
    JCalc(string name){
        Init(name);
    }

    CrgUnitType * GetCrgUnitTypeByName(string);

    void Init(string );

    vector <double> GetJ (CrgUnit & one, CrgUnit & two);
    CrgUnit DefineCrgUnit(vec pos, matrix orient, string name);
private:

    struct JCalcData{
            /// the CrgUnitType of the first (will contain coordinates of the still molecules
        CrgUnitType *_type1;
        CrgUnitType *_type2;

        ///variables copied from _type1 and _type2
        pair <vector <unsigned int>, vector <unsigned int> > _orblabels;

        /// the moved coordinates
        mol_and_orb _mol1;
        mol_and_orb _mol2;
        /// the moved orbitals
        orb         _orb1;
        orb         _orb2;
        /// the fock matrix herself
        fock        _fock;

        basis_set   _indo;
    };
    
    /// a list of charge transport units
    vector <CrgUnitType *> _listCrgUnitType;
    /// a map of charge unit type names to their index in the list
    map <string, CrgUnitType *> _mapCrgUnitByName;
    /// map the fock calculators to the pair of integers representing the corresponding charge unit types
    map <pair<CrgUnitType *, CrgUnitType *> , JCalcData *> _maplistfock;

    /// Enter data for Crg Unit Type
    void ParseCrgUnitType(xmlDocPtr doc, xmlNodePtr cur );
    /// initialise a JCAlcDAta type
    JCalcData * InitJCalcData(CrgUnitType *, CrgUnitType * );
};

#endif	/* _JCALC_H */

