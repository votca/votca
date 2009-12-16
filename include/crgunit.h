/* 
 * File:   crgunit.h
 * Author: james
 *
 * Created on 27 June 2007, 09:04
 */

#ifndef _CRGUNIT_H
#define	_CRGUNIT_H

#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include "global.h"
#include "crgunittype.h"

using namespace std;

/**
    \brief information about a CrgUnit
 
    The CrgUnit class describes a charge unit. It stores information like the id, the position,
 * orientation, type of the charge transport unit. It 

 */
class CrgUnit {
public:

    CrgUnit() {
    };

    ~CrgUnit();

    CrgUnit(vector <vec> positions, vector <vec> norms, vector <vec> planes,
            const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId);
    CrgUnit(vec pos, vec plane1, vec norm1, CrgUnitType * type);
    
    /*CrgUnit(CrgUnit & acrg){
        copyCrgUnit(acrg);
    }*/
    
    void copyCrgUnit(CrgUnit & acrg);
    void copyCrgUnit(CrgUnit & acrg, const int & id);
    const double& GetNRG() const;

    void setNRG(const double& a);

    const unsigned int& GetId() const;
    void SetId(const int& i);

    const unsigned int& GetTypeID() const;
    CrgUnitType* GetType() const;

    const unsigned int& GetMolId() const;

    int GetN();

    const vector <vec>::iterator GetPosFirst();
    const vector <vec>::iterator GetPosLast();

    void SetCom(vec& com);
    void SetPos(const int& i, const vec& pos);
    void SetNorm(const int& i, const vec& pos);
    void SetPlane(const int& i, const  vec& pos);

    mol_and_orb * rotate_translate_beads();

    // this is the function called from the rate calculator on already initialised molecules
    void rot_two_mol(CrgUnit & two, mol_and_orb & mol1, mol_and_orb & mol2);

    vec GetPos(const int & i);
    vec GetNorm(const int & i);
    vec GetPlane(const int & i);
    const vec & GetCom() const;

    void shift(const vec & displ);

    void rotate(matrix mat);

private:
    /// the centre of mass
    vec _com;
    /// the orientation
    //   matrix _orient ;
    /// the ID 
    unsigned int _id;
    /// the type
    string _name;
    /// the molecule index
    unsigned int _molid;
    /// also the type, but as an int
    //unsigned int _typeid;
    ///vector of coms of monomers
    vector < vec > _positions;
    ///vector of coms of monomers
    vector < vec > _norms;
    ///vector of coms of monomers
    vector < vec > _planes;
    //vector of coordinates (if on a ghost cell)
    //vector < vector < vec > > _altpos;
    // a vector of the alternative centre of masses
    //vector < vec > _altcom;
    /// a reference to the crgunittype
    CrgUnitType * _type;
    /// the energy at this site
    double _energy;

    vector <vec>::iterator GetPositions() {
        return _positions.begin();
    }

    vector <vec>::iterator GetNorms() {
        return _norms.begin();
    }

    vector <vec>::iterator GetPlanes() {
        return _planes.begin();
    }

    vector <vec> shift_pos(const vec & a);
};

inline void CrgUnit::copyCrgUnit(CrgUnit & acrg, const int & id) {
    copyCrgUnit(acrg);
    _id = id;
}

inline const double& CrgUnit::GetNRG() const {
    return _energy;
}

inline void CrgUnit::setNRG(const double& a) {
    _energy = a;
}

inline const unsigned int& CrgUnit::GetId() const {
    return _id;
}

inline void CrgUnit::SetId(const int& i) {
    _id = i;
}

inline const unsigned int& CrgUnit::GetTypeID() const {
    return _type->GetId();
}

inline CrgUnitType* CrgUnit::GetType() const {
    return _type;
}

inline const unsigned int& CrgUnit::GetMolId() const {
    return _molid;
}

inline int CrgUnit::GetN() {
    if (_norms.size() != _planes.size() || _norms.size() != _positions.size()) {
        cerr << "ERror in the crg unit of id:" << _id << endl;
    }
    return _norms.size();
}

inline const vector <vec>::iterator CrgUnit::GetPosFirst() {
    return _positions.begin();
}

inline const vector <vec>::iterator CrgUnit::GetPosLast() {
    return _positions.end();
}

inline void CrgUnit::SetCom(vec& com) {
    _com = com;
};

inline void CrgUnit::SetPos(const int& i, const vec& pos) {
    if ( i >= _positions.size()){
        _positions.resize(i+1);
        _norms.resize(i+1);
        _planes.resize(i+1);
    }
    _positions[i] = pos;
};
inline void CrgUnit::SetNorm(const int& i, const vec& pos) {
    if ( i >= norms.size()){
        _positions.resize(i+1);
        _norms.resize(i+1);
        _planes.resize(i+1);
    }
    _norms[i] = pos;
};
inline void CrgUnit::SetPlane(const int& i, const vec& pos) {
    if ( i >= _planes.size()){
        _positions.resize(i+1);
        _norms.resize(i+1);
        _planes.resize(i+1);
    }
    _planes[i] = pos;
};

inline vec CrgUnit::GetPos(const int & i) {
    return _positions[i];
}

inline vec CrgUnit::GetNorm(const int & i) {
    return _norms[i];
}

inline vec CrgUnit::GetPlane(const int & i) {
    return _planes[i];
}

inline const vec & CrgUnit::GetCom() const {
    return _com;
}

inline void CrgUnit::shift(const vec & displ) {
    vector < vec> ::iterator it_vec;
    for (it_vec = _positions.begin(); it_vec != _positions.end(); ++it_vec) *it_vec = *it_vec + displ;
    _com += displ;
}

inline vector <vec> CrgUnit::shift_pos(const vec & a) {
    vector <vec>::iterator it_pos;
    vector <vec> res;
    for (it_pos = _positions.begin(); it_pos != _positions.end(); ++it_pos) {
        vec b = (*it_pos) + a;
        res.push_back(b);
    }
    return res;
}

#endif	/* _CRGUNIT_H */

