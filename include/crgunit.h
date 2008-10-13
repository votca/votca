/* 
 * File:   crgunit.h
 * Author: james
 *
 * Created on 27 June 2007, 09:04
 */

#ifndef _CRGUNIT_H
#define	_CRGUNIT_H

#include <tools/vec.h>
#include <tools/matrix.h>
#include "global.h"
#include "crgunittype.h"

using namespace std;

/**
    \brief information about a CrgUnit
 
    The CrgUnit class describes a charge unit. It stores information like the id, the position,
 * orientation, type of the charge transport unit. It 

*/
class CrgUnit{
public:
    CrgUnit(){};
    
    ~CrgUnit(){
        
        
        cout << "[crgunit.h]: Calling the crgunit destructor." << endl;

        _planes.clear();
        _positions.clear();
        _norms.clear();
    }
    
    CrgUnit(vector <vec> positions, vector <vec> norms, vector <vec> planes, 
            const unsigned int & id,  CrgUnitType * type , const unsigned int & molId){
        vector <vec>::iterator it_pos;
        vector <vec>::iterator it_norm;
        vector <vec>::iterator it_planes;
        _com = vec(0.,0.,0.);
        double toll = 1E-6;
        vec xvec (1.0,0.,0.);
        vec yvec (0.,1.,0.);
        for (it_pos = positions.begin() , it_norm = norms.begin(), it_planes = planes.begin()
        ; it_pos != positions.end(); ++it_pos, ++it_norm , ++it_planes){
            _positions.push_back(*it_pos);
            if (abs(*it_norm) > toll  && abs(*it_planes ) > toll){
                _norms.push_back(*it_norm);
                _planes.push_back(*it_planes);
            }
            else {
                _norms.push_back(xvec);
                _planes.push_back(yvec);
            }
            _com = _com + *it_pos;
            
        }
        _com /= positions.size();
        _molid = molId;
        _id  = id;
        _type = type;        
        
        // initialise variables derived from CrgUnitType
        _energy = _type -> GetEnergy();
    }
    
    void copyCrgUnit( CrgUnit & acrg, const int & id){
        vector <vec>::iterator it_pos;
        vector <vec>::iterator it_norm;
        vector <vec>::iterator it_planes;
        _com = acrg._com;
        _positions.clear();
        _norms.clear();
        _planes.clear();
        for (it_pos = (acrg._positions).begin() , 
             it_norm = (acrg._norms).begin(), 
             it_planes = (acrg._planes).begin()
        ; it_pos != (acrg._positions).end(); ++it_pos, ++it_norm , ++it_planes){
            _positions.push_back(*it_pos);
            _norms.push_back(*it_norm);
            _planes.push_back(*it_planes);
            
        }
        _molid = acrg._molid;
        _id  = id;
        _type = acrg._type;        
        
        // initialise variables derived from CrgUnitType
        _energy = acrg._energy;
    
    }
    
    const double& GetNRG() const{
        return _energy;
    }
    
    void setNRG(const double& a){
        _energy=a;
    }
    const unsigned int& GetId() const{
        return _id;
    }
    
    void SetId (const int& i) {
        _id = i;
    }
    
    const unsigned int& GetTypeID() const{
        return _type->GetId();
    }
    
    CrgUnitType* GetType() const{
        return _type;
    }
    
    const unsigned int& GetMolId() const{
        return _molid;
    }
    
    int GetN(){
        if (_norms.size() != _planes.size() || _norms.size() != _positions.size()){
            cerr << "ERror in the crg unit of id:" << _id <<endl;
        }
        return _norms.size();
    }
    
    const vector <vec>::iterator  GetPosFirst () {
        return _positions.begin();
    }
    const vector <vec>::iterator  GetPosLast () {
        return _positions.end();
    }
    void SetCom(vec& com){
        _com = com;
    };
    void SetPos(const int& i, vec& pos){
        _positions[i] = pos;
    };

    mol_and_orb * rotate_translate_beads( ){
        mol_and_orb *a = new mol_and_orb;
        basis_set *_indo = new basis_set;
        orb * _orb = new orb;
        a->define_bs(*_indo);
        a->cp_atompos(_type->GetCrgUnit() );
        a->cp_atoms  (_type->GetCrgUnit() );
        _orb->init_orbitals_stripped(_type->GetOrb());
        a->assign_orb(_orb);
        
        _type->rotate_each_bead( this->GetPositions() , this->GetNorms(), 
                 this->GetPlanes(), a);
        return a;
    }
    
    ///strips the extra copies of positions correspondinig to ghost atoms
   /* void strip_bc(){
        clearListList(_altpos);
        _altcom.clear();
    }*/
    
    /// generates extra copies of the positions corresponding to "ghost" atoms
   /* void apply_bc ( const vec & bc, const double & d){
        strip_bc();
        
        if ( _com.getX() > (bc.getX()-d) ){
            vec displ (-bc.getX(), 0. ,0. );
            vector <vec> pos = shift_pos(displ);
            _altpos.push_back(pos);
            vec com2 = _com + displ;
            _altcom.push_back(com2);
        }
        else if ( _com.getX() < d){
            vec displ ( bc.getX(), 0. ,0. );
            vector <vec> pos = shift_pos(displ);
            _altpos.push_back(pos);
            vec com2 = _com + displ;
            _altcom.push_back(com2);
        }
        if ( _com.getY() > (bc.getY()-d) ){
            vec displ (0. , -bc.getY(), 0. );
            vector <vec> pos = shift_pos(displ);
            _altpos.push_back(pos);
            vec com2 = _com + displ;
            _altcom.push_back(com2);
        }
        else if ( _com.getY() < d){
            vec displ (0. , bc.getY(), 0. );
            vector <vec> pos = shift_pos(displ);
            _altpos.push_back(pos); 
            vec com2 = _com + displ;
            _altcom.push_back(com2);      
        }
        if ( _com.getZ() > (bc.getZ()-d) ){
            vec displ (0. , 0, -bc.getZ() );
            vector <vec> pos = shift_pos(displ);
            _altpos.push_back(pos);
            vec com2 = _com + displ;
            _altcom.push_back(com2);
        }
        else if ( _com.getZ() < d){
            vec displ (0. , 0, bc.getZ() );
            vector <vec> pos = shift_pos(displ);
            _altpos.push_back(pos); 
            vec com2 = _com + displ;
            _altcom.push_back(com2);       
        }
        
    }*/

    ///find the distance between two charge units
    /*double d_min(CrgUnit &a){
        return abs( displ_min(a));
    }*/
        
    // this is the function called from the rate calculator on already initialised molecules
    void rot_two_mol(CrgUnit & two, mol_and_orb & mol1, mol_and_orb & mol2){
        
        /*vector <vec> ::iterator  com1;
        vector <vec> ::iterator  com2;
        vector <vector < vec > > :: iterator pos1;
        vector <vector < vec > >:: iterator pos2;*/
        
        vector <vec>::iterator bestpos1;
        vector <vec>::iterator bestpos2;
        
        bestpos1 = _positions.begin();
        bestpos2 = two._positions.begin();
        
        /*double d = abs( two._com - _com );
        
        for ( com2 = (two._altcom).begin(), pos2 = (two._altpos).begin() ;
        com2 != (two._altcom).end() ; ++com2, ++pos2 ){
            vec displ = *com2 - _com;
            double d_t = abs( displ);
            if ( d_t < d ) {
                bestpos1 = _positions.begin();
                bestpos2 = pos2 -> begin();
                d=d_t;
            }
        }
        
        for ( com1 = (_altcom).begin(), pos1 = (_altpos).begin() ;
        com1 != (_altcom).end() ; ++com1, ++pos1 ){
            vec displ = two._com - *com1;
            double d_t = abs( displ);
            if ( d_t < d ) {
                bestpos1 = pos1 -> begin();
                bestpos2 = two._positions.begin();
                d=d_t;
            }
        }
        
        for ( com1 = _altcom.begin(), pos1 = _altpos.begin() ; 
        com1 != _altcom.end() ; ++com1, ++pos1 ){
            for ( com2 = (two._altcom).begin(), pos2 = (two._altpos).begin() ;
            com2 != (two._altcom).end() ; ++com2, ++pos2 ){
                vec displ = *com2 - *com1;
                double d_t = abs( displ);
                if ( d_t < d ) {
                    bestpos1 = pos1 -> begin();
                    bestpos2 = pos2 -> begin();
                    d=d_t;
                }
            }
        }*/
        
        _type->rotate_each_bead(bestpos1,this->GetNorms(), this->GetPlanes(), &mol1 );
        (two._type)->rotate_each_bead(bestpos2, two.GetNorms(), two.GetPlanes(), &mol2);
    }
    
    vec GetPos(const int & i){
        return _positions[i];
    }
    vec GetNorm(const int & i){
        return _norms[i];
    }
    vec GetPlane(const int & i){
        return _planes[i];
    }
    const vec & GetCom () const{
        return _com;
    }
    vec GetComInBox (const vec & bc) const{
        vec com=GetCom();
        while (com.getX() > bc.getX() ) com.setX(com.getX() - bc.getX() );
        while (com.getX() < 0 )    com.setX(com.getX() + bc.getX() );
        while (com.getY() > bc.getY() ) com.setY(com.getY() - bc.getY() );
        while (com.getY() < 0 )    com.setY(com.getY() + bc.getY() );
        while (com.getZ() > bc.getZ() ) com.setZ(com.getZ() - bc.getZ() );
        while (com.getZ() < 0 )    com.setZ(com.getZ() + bc.getZ() );
        return com;
    }
    
    void shift(const vec & displ) {
        _com = _com + displ;
        vector < vec> ::iterator it_vec;
        //for (it_vec = _altcom.begin() ; it_vec != _altcom.end(); ++it_vec) *it_vec = *it_vec + displ;
        for (it_vec = _positions.begin() ; it_vec != _positions.end(); ++it_vec) *it_vec = *it_vec + displ;
        /*vector < vector < vec > > :: iterator it_altpos;
        for (it_altpos = _altpos.begin() ; it_altpos != _altpos.end() ; ++it_altpos){
            for (it_vec =  it_altpos->begin(); it_vec != it_altpos->end() ; ++it_vec ){
                *it_vec = *it_vec + displ;
            }
        }*/
    }
    
    void rotate(matrix mat){
        vector <vec> ::iterator it_norm = _norms.begin();
        vector <vec> ::iterator it_plane = _planes.begin();
        for ( ; it_norm != _norms.end() ; ++it_norm, ++it_plane){
            vec a = *it_norm;
            vec b = *it_plane;
            *it_norm  = mat * a;
            *it_plane = mat * b;
        }
        
        vec a = mat * _com ;
        _com = a;
        vector < vec> ::iterator it_vec;
        /*for (it_vec = _altcom.begin() ; it_vec != _altcom.end(); ++it_vec) {
            a = *it_vec;
            *it_vec = mat * (*it_vec);
        }*/
        for (it_vec = _positions.begin() ; it_vec != _positions.end(); ++it_vec) {
            a = *it_vec;
            *it_vec = mat * (*it_vec);
        }
        /*vector < vector < vec > > :: iterator it_altpos;
        for (it_altpos = _altpos.begin() ; it_altpos != _altpos.end() ; ++it_altpos){
            for (it_vec =  it_altpos->begin(); it_vec != it_altpos->end() ; ++it_vec ){
                a = *it_vec;
                *it_vec = mat * (*it_vec);
            }
        }*/
    }
    
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
    
    vector <vec>::iterator  GetPositions () {
        return _positions.begin();
    }
    vector <vec>::iterator  GetNorms () {
        return _norms.begin();
    }
    vector <vec>::iterator GetPlanes () {
        return _planes.begin();
    }
    
    
    vector <vec> shift_pos (const vec & a){
        vector <vec>::iterator it_pos; 
        vector <vec> res;
        for (it_pos = _positions.begin(); it_pos != _positions.end() ;++it_pos){
            vec b = (*it_pos) + a;
            res.push_back(b);
        }
        return res;
    }
} ;

#endif	/* _CRGUNIT_H */

