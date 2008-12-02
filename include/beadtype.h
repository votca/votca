/* 
 * File:   beadtype.h
 * Author: ruehle
 *
 * Created on November 26, 2008, 12:44 PM
 */

#ifndef _BEADTYPE_H
#define	_BEADTYPE_H

#include <string>
using namespace std;

class  BeadType {
public:
    BeadType(int id, const string &name)
    : _id(id), _name(name) {}
    
    const int &getId() const { return _id; }
    const string &getName() const { return _name; }
    
private:
    int _id;
    string _name;
};

#endif	/* _BEADTYPE_H */

