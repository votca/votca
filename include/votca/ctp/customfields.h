/* 
 * File:   customfields.h
 * Author: ruehle
 *
 * Created on July 21, 2011, 11:43 AM
 */

#ifndef __VOTCA_MD2QM_CUSTOMFIELDS_H
#define	__VOTCA_MD2QM_CUSTOMFIELDS_H

#include <map>
#include <string>

class CustomFields
{
public:
    void setDouble(std::string key, double value) { _double_values[key] = value; }
    double getDouble(std::string key) { return _double_values[key]; }
    std::map<std::string, double> &DoubleValues() { return _double_values; }
    bool DoubleExists(std::string key) {
   return _double_values.find(key) != _double_values.end();
} 
protected:
    std::map<std::string, double> _double_values;
};


#endif	/* CUSTOMFIELDS_H */

