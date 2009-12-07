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


#include "crgunittype.h"
#include "votca/tools/tokenizer.h"

class JCalc{

public:

    JCalc(){

    }
    JCalc(string name){
        Init(name);
    }

    void Init(string );
private:
    vector <CrgUnitType *> _listCrgUnitType;



void ParseCrgUnitType(xmlDocPtr doc, xmlNodePtr cur );

};

#endif	/* _JCALC_H */

