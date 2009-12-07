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

#include "crgunittype.h"
#include "votca/tools/tokenizer.h"

using namespace boost;

class JCalc{

public:

    JCalc(){

    }
    JCalc(string name){
        Init(name);
    }

    void Init(string );
private:
/// a list of charge transport units
    vector <CrgUnitType *> _listCrgUnitType;
    /// a map of charge unit type names to their index in the list
    map <string, CrgUnitType *> _mapCrgUnitByName;



void ParseCrgUnitType(xmlDocPtr doc, xmlNodePtr cur );

};

#endif	/* _JCALC_H */

