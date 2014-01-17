#ifndef _TOP2Map__H
#define _TOP2Map__H


#include <votca/ctp/topology.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/logger.h>
#include <boost/algorithm/string.hpp>
#include <votca/tools/vec.h>
#include <boost/format.hpp>


namespace votca { namespace ctp {
    namespace ba = boost::algorithm;
    using namespace std;
    
class TOP2Map : public QMTool
{
public:

    TOP2Map() { };
   ~TOP2Map() { };
   
    string Identify() { return "top2map"; }
    // run sequence
    void   Initialize(Property *options);
    bool   Evaluate();

    // helpful guys
    void readPDB();
    void readGRO();
    
    // make top file
    void top2txt();
    
    // error formating function
    void error1(string line){ cout << endl; throw runtime_error(line); }
    
private:
    string      _input_pdb;
    string      _input_gro;
    
    bool        _has_pdb;
    bool        _has_gro;
    
    int         _numMol;
    
    Topology    _top;
};

void TOP2Map::Initialize(Property* options) {
    // read options    
    string key = "options.top2map.";

    // set boolean constants to false
    _has_pdb = false;
    _has_gro = false;
    
    // find PDB, then GRO, then error
    if ( options->exists(key+"pdb") ){
        _input_pdb      = options->get(key+"pdb").as<string> ();
        _has_pdb        = true;
        cout << endl << "... ... PDB file: \t" << _input_pdb;
    }
    else if ( options->exists(key+"gro") ){
        _input_gro      = options->get(key+"gro").as<string> ();
        _has_gro        = true;
        cout << endl << "... ... GRO file: \t" << _input_gro;
    }
    else{
        error1( "... ... Error. Unsupported input file format. \n"
                "... ... Supported formats: pdb, gro \n");     
    }
    
    if ( options->exists(key+"mol_num") ){
        _numMol = options->get(key+"mol_num").as<int> ();
        cout << "\n" "... ... User defined num of mols: \t" << _numMol;
    }
    else{_numMol = 1; cout <<"\n" "... ... Default num of mols: 1 \n";}
}

void TOP2Map::top2txt(){
    Topology * _topPtr = &_top;
//    Segment * _seg = _topPtr->Molecules().back()->Segments().back();

}

bool TOP2Map::Evaluate() {
    top2txt();

}

void TOP2Map::readPDB(){
   cout << endl << "... ... Assuming: PDB";
    
    // set molecule >> segment >> fragment
    // reconnect them all
    Topology * _topPtr = 0;
    _topPtr = &_top;
    
    Molecule * _molPtr = 0;
    // direct
    _molPtr = _topPtr->AddMolecule("M1");
                // inverse
                _molPtr->setTopology(_topPtr);
    
    Segment  * _segPtr  = 0;
    // direct
    _segPtr = _topPtr->AddSegment("S1");
               _molPtr->AddSegment(_segPtr);
               // inverse
                _segPtr->setTopology(_topPtr);
                _segPtr->setMolecule(_molPtr);

    // try: read PDB file
    std::ifstream _file( _input_pdb.c_str());
    if (_file.fail()) {
        error1( "... ... Can not open: " + _input_pdb + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    }
    else{
        cout << endl << 
                ("... ... File " + _input_pdb + ""
                 " was opened successfully.\n");
    }

    // read PDB line by line
    string _line;
    
    // counters for loops
//    int _atom_id = 0;
    int _newResNum = 0;

    while ( std::getline(_file, _line,'\n') ){
        if(     boost::find_first(_line, "ATOM"  )   || 
                boost::find_first(_line, "HETATM")      
                ){
            
            //      according to PDB format
            string _recType    (_line,( 1-1),6); // str,  "ATOM", "HETATM"
            string _atNum      (_line,( 7-1),6); // int,  Atom serial number
            string _atName     (_line,(13-1),4); // str,  Atom name
            string _atAltLoc   (_line,(17-1),1); // char, Alternate location indicator
            string _resName    (_line,(18-1),4); // str,  Residue name
            string _chainID    (_line,(22-1),1); // char, Chain identifier
            string _resNum     (_line,(23-1),4); // int,  Residue sequence number
            string _atICode    (_line,(27-1),1); // char, Code for insertion of res
            string _x          (_line,(31-1),8); // float 8.3 ,x
            string _y          (_line,(39-1),8); // float 8.3 ,y
            string _z          (_line,(47-1),8); // float 8.3 ,z
            string _atOccup    (_line,(55-1),6); // float  6.2, Occupancy
            string _atTFactor  (_line,(61-1),6); // float  6.2, Temperature factor
            string _segID      (_line,(73-1),4); // str, Segment identifier
            string _atElement  (_line,(77-1),2); // str, Element symbol
            string _atCharge   (_line,(79-1),2); // str, Charge on the atom

            ba::trim(_recType);
            ba::trim(_atNum);
            ba::trim(_atName);
            ba::trim(_atAltLoc);
            ba::trim(_resName);
            ba::trim(_chainID);
            ba::trim(_resNum);
            ba::trim(_atICode);
            ba::trim(_x);
            ba::trim(_y);
            ba::trim(_z);
            ba::trim(_atOccup);
            ba::trim(_atTFactor);
            ba::trim(_segID);
            ba::trim(_atElement);
            ba::trim(_atCharge);
            
            double _xd(0),_yd(0),_zd(0);
            int _resNumInt(0); 
            
            try
            {
            _xd = boost::lexical_cast<double>(_x);
            _yd = boost::lexical_cast<double>(_y);
            _zd = boost::lexical_cast<double>(_z);
            _resNumInt = boost::lexical_cast<int>(_resNum);
            }
            catch(boost::bad_lexical_cast &)
            {
                error1( "... ... Can not convert PDB coord line!\n"
                        "... ... Atom number: " + _atNum + "\n"
                        "... ... Make sure this line is PDB style\n");
            }
            
            vec r(_xd , _yd , _zd);

            // set fragment
            // reconnect to topology, molecule, segment
            Fragment * _fragPtr = 0;
            // make new frag for new res number
            // otherwise use last created
            if ( _newResNum != _resNumInt ){

                _newResNum = _resNumInt;
                string _newResName = _resName+'_'+_resNum;
                
                // direct
                _fragPtr = _topPtr->AddFragment(_newResName);
                           _molPtr->AddFragment(_fragPtr);
                           _segPtr->AddFragment(_fragPtr);
                          // inverse
                          _fragPtr->setTopology(_topPtr);
                          _fragPtr->setMolecule(_molPtr);
                          _fragPtr->setSegment(_segPtr);        
            }
            else{
                _fragPtr = _topPtr->Fragments().back();
            }
            if (_fragPtr==0) {error1("Zero pointer in GRO reader. Why?");}
                        
            // set atom
            // reconnect to topology, molecule, segment, fragment
            Atom * _atmPtr = 0;
            // direct
            _atmPtr = _topPtr->AddAtom(_atName);
                      _molPtr->AddAtom(_atmPtr);
                      _segPtr->AddAtom(_atmPtr);
                     _fragPtr->AddAtom(_atmPtr);
                      // inverse
                      _atmPtr->setTopology(_topPtr);
                      _atmPtr->setMolecule(_molPtr);        
                      _atmPtr->setSegment(_segPtr);
                      _atmPtr->setFragment(_fragPtr);
                      
            _atmPtr->setResnr        (_resNumInt);
            _atmPtr->setResname      (_resName);
            _atmPtr->setPos          (r);
        }
    }
    
    return;
}

void TOP2Map::readGRO(){
    cout << endl << "... ... Assuming: GRO";

    // set molecule >> segment >> fragment
    // reconnect them all
    Topology * _topPtr = 0;
    _topPtr = &_top;
    
    Molecule * _molPtr = 0;
    // direct
    _molPtr = _topPtr->AddMolecule("M1");
                // inverse
                _molPtr->setTopology(_topPtr);
    
    Segment  * _segPtr  = 0;
    // direct
    _segPtr = _topPtr->AddSegment("S1");
               _molPtr->AddSegment(_segPtr);
               // inverse
                _segPtr->setTopology(_topPtr);
                _segPtr->setMolecule(_molPtr);

    // try: read GRO file
    std::ifstream _file( _input_gro.c_str());
    if (_file.fail()) {
        error1( "... ... Can not open: " + _input_gro + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    }
    else{
        cout << endl << 
                ("... ... File " + _input_gro + ""
                 " was opened successfully.\n");
    }

    // read GRO line by line
    string _line;
    
    // counters for loops
    int _newResNum = -1; // res reference
    int _atTotl = 0;  // total num of atoms in GRO
    int _atCntr = 0;  // atom number counter
    
    // GRO: first two lines are tech specs -> ignore them
    // ignore first line, it's a comment
    std::getline(_file, _line,'\n');

    // GRO check: if second line can cast to int, then ok

    try
    {   
        // first line, number of atoms in XYZ
        std::getline(_file, _line,'\n');
        ba::trim(_line);
        _atTotl = boost::lexical_cast<int>(_line);
    }
    catch(boost::bad_lexical_cast &)
    {
        error1( "... ... Bad GRO file format!\n"
                "... ... First two line must contain technical specs.\n");
    }

    // actual loop
    while ( std::getline(_file, _line,'\n') ){
        if (_atCntr < _atTotl){
            
            string _resNum     (_line, 0,5); // int,  Residue number
            string _resName    (_line, 5,5); // str,  Residue name
            string _atName     (_line,10,5); // str,  Atom name
            string _atNum      (_line,15,5); // int,  Atom number
            string _x          (_line,20,8); // float 8.3 ,x
            string _y          (_line,28,8); // float 8.3 ,y
            string _z          (_line,36,8); // float 8.3 ,z
            
            ba::trim(_atNum);
            ba::trim(_atName);
            ba::trim(_resNum);
            ba::trim(_resName);
            ba::trim(_x);
            ba::trim(_y);
            ba::trim(_z);
            
            // try cast
            int _resNumInt(0),_atNumInt(0);
            double _xd(0),_yd(0),_zd(0);
            try
            {
                _resNumInt = boost::lexical_cast<int>(_resNum);
                _atNumInt  = boost::lexical_cast<int>(_atNum);

                _xd = boost::lexical_cast<double>(_x);
                _yd = boost::lexical_cast<double>(_y);
                _zd = boost::lexical_cast<double>(_z);
            }
            catch (boost::bad_lexical_cast &)
            {
                error1( "... ... Can not convert GRO coord line!\n"
                        "... ... Atom number: " + _atNum + "\n"
                        "... ... Make sure this line is GRO style\n");
            }
            
            vec r(_xd , _yd , _zd);
                
            // set fragment
            // reconnect to topology, molecule, segment
            Fragment * _fragPtr = 0;
            // make new frag for new res number
            // otherwise use last created
            if ( _newResNum != _resNumInt ){

                _newResNum = _resNumInt;
                string _newResName = _resName+'_'+_resNum;
                
                // direct
                _fragPtr = _topPtr->AddFragment(_newResName);
                           _molPtr->AddFragment(_fragPtr);
                           _segPtr->AddFragment(_fragPtr);
                          // inverse
                          _fragPtr->setTopology(_topPtr);
                          _fragPtr->setMolecule(_molPtr);
                          _fragPtr->setSegment(_segPtr);        
            }
            else{
                _fragPtr = _topPtr->Fragments().back();
            }
            if (_fragPtr==0) {error1("Zero pointer in GRO reader. Why?");}
                        
            // set atom
            // reconnect to topology, molecule, segment, fragment
            Atom * _atmPtr = 0;
            // direct
            _atmPtr = _topPtr->AddAtom(_atName);
                      _molPtr->AddAtom(_atmPtr);
                      _segPtr->AddAtom(_atmPtr);
                     _fragPtr->AddAtom(_atmPtr);
                      // inverse
                      _atmPtr->setTopology(_topPtr);
                      _atmPtr->setMolecule(_molPtr);        
                      _atmPtr->setSegment(_segPtr);
                      _atmPtr->setFragment(_fragPtr);
        
            _atmPtr->setResnr        (_resNumInt);
            _atmPtr->setResname      (_resName);
            _atmPtr->setPos          (r);
        
        }
        _atCntr++;
    }
    
    return;
}

}}


#endif
