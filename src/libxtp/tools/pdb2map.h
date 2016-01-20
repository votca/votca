#ifndef _PDB2Map_H
#define _PDB2Map_H


#include <votca/xtp/topology.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/logger.h>
#include <boost/algorithm/string.hpp>
#include <votca/tools/vec.h>
#include <boost/format.hpp>


namespace votca { namespace xtp {
    namespace ba = boost::algorithm;
    using namespace std;
    
class PDB2Map : public QMTool
{
public:

    PDB2Map() { };
   ~PDB2Map() { };
   
    string Identify() { return "pdb2map"; }
    
    // read options    
    void   Initialize(Property *options);
    
    // make xml
    bool   Evaluate();
    
    // helpful guys
    void readPDB();
    void readGRO();
    void readXYZ();
    void setTopologies();
    void adaptQM2MD();
    void topMdQm2xml();
    
    void error1(string line){ cout << endl; throw runtime_error(line); };
 
    
private:
    string      _input_pdb;
    string      _input_gro;
    string      _input_xyz;
    
    string      _output_xml;
    
    bool        _has_pdb;
    bool        _has_gro;
    bool        _has_xyz;    

    bool        _has_md;    
    bool        _has_qm;    
    
    bool        _can_convert_md2qm;
    bool        _QM2MDcompatible;
    
    Topology    _MDtop;
    Topology    _QMtop;

   // element:mass map
    map <string,int> el2mass;

};

void PDB2Map::Initialize(Property* options) 
{   

    // update options with the VOTCASHARE defaults   
    //    UpdateWithDefaults( options );
    
    // fill in periodic table
    el2mass["H"]        = 1;
    el2mass["B"]        = 10;
    el2mass["C"]        = 12;
    el2mass["N"]        = 14;
    el2mass["O"]        = 16;
    el2mass["F"]        = 19;
    el2mass["Al"]       = 27;
    el2mass["Si"]       = 28;
    el2mass["P"]        = 31;
    el2mass["S"]        = 32;
    el2mass["Cl"]       = 35;
    el2mass["Ir"]       = 192;
    
    // file exists
    _has_pdb = false;
    _has_gro = false;
    _has_xyz = false;
    
    // top exists
    _has_qm  = false;
    _has_md  = false;
    
    // can convert from MD to QM, if QM is not found
    _can_convert_md2qm = false;
    
    // read options    
    string key = "options.pdb2map.";
    
    // find PDB, then GRO, then error
    if ( options->exists(key+"pdb") ){
        _input_pdb      = options->get(key+"pdb").as<string> ();
        _has_pdb        = true;
        cout << endl << _input_pdb;
        cout << endl << "... ... PDB input specified: \t" << _input_pdb;
    }
    else if ( options->exists(key+"gro") ){
        _input_gro      = options->get(key+"gro").as<string> ();
        _has_gro        = true;
        cout << endl << "... ... GRO input specified: \t" << _input_gro;
    }
    else{
        error1("... ... No MD(PDB,GRO) file provided. "
                                "\n... ... Tags: pdb, gro");     
    }
    
    // find XYZ, then error
    if ( options->exists(key+"xyz") ){
        _input_xyz      = options->get(key+"xyz").as<string> ();
        _has_xyz        = true;
        cout << endl << "... ... XYZ input specified: \t" << _input_xyz;
    }
    else if (_has_pdb){
            cout << endl << "... ... *** No QM(XYZ) file provided. Tags: xyz\n"
                  "... ... BUT I can make a map from PDB only. Continue.";
    }
    else{
        error1("... ... No QM(XYZ) file provided. Tags: xyz "
                                "\n... ... No QM(PDB) as substitute.");    
    }

    // find XML or generate it
    if ( options->exists(key+"xml") ){
        _output_xml = options->get(key+"xml").as<string> ();
        cout << endl << "... ... XML output specified: \t" << _output_xml;
        }
        else{
                _output_xml = "system_output.xml";
                cout << endl << "... ... *** No XML output specified.";
                cout << endl << "... ... Default XML is: \t" << _output_xml;
        }

}

bool PDB2Map::Evaluate() {
    
    setTopologies();
    
    // if has XYZ
    // XYZ has no residues
    // adaptQM2MD() maps MD residues to QM top
    // in order to iterate properly
    if (_has_xyz){adaptQM2MD();}

    topMdQm2xml();
    
//    LOG( logINFO, _log ) << "Reading from: " << _input_file << flush;    
//    cout << _log;
    return true;
}

void PDB2Map::setTopologies(){
    if (_has_pdb){
        readPDB();
        _has_md = true;
    }
    else if (_has_gro){
        readGRO();
        _has_md = true;
    }
    else{
        error1("... ... Bad MD topology.\n"
               "... ... Please supply PDB or GRO file");
    }
    
    if (_has_xyz){
        readXYZ();
        _has_qm = true;
    }
    else if (_can_convert_md2qm){
        _QMtop = _MDtop;
        _has_qm = true;
    }
    else{
        error1("... ... Bad MD topology.\n"
               "... ... Please supply XYZ file or PDB with chemical elements");
    }
}


void PDB2Map::adaptQM2MD(){
    
    // check if atom number is the same in QM and MD
    int numMDatoms = _MDtop.getMolecule(1)->NumberOfAtoms();
    int numQMatoms = _QMtop.getMolecule(1)->NumberOfAtoms();
    
    _QM2MDcompatible = (numMDatoms == numQMatoms) ? true : false;
    
    // if so, proceed
    if (_QM2MDcompatible){
        Molecule * MDmolecule = _MDtop.getMolecule(1);
        Molecule * QMmolecule = _QMtop.getMolecule(1);
        
        vector<Segment*> MDsegments = MDmolecule->Segments();
        vector<Segment*> QMsegments = QMmolecule->Segments();
        
        vector<Segment*>::iterator MDSegIt;
        vector<Segment*>::iterator QMSegIt;
        
        for(MDSegIt = MDsegments.begin(), QMSegIt = QMsegments.begin();
                MDSegIt < MDsegments.end();
                MDSegIt++, QMSegIt++  )
        {
            Fragment * QMfragment = 0;
            
            vector<Atom*> MDSegAtoms = (*MDSegIt)->Atoms();
            vector<Atom*> QMSegAtoms = (*QMSegIt)->Atoms();
            
            vector<Atom*>::iterator MDSegAtIt;
            vector<Atom*>::iterator QMSegAtIt;
            
            int old_res_num(-1),new_res_num(-1);
            string res_name = "bad_wolf";
            
            for(MDSegAtIt = MDSegAtoms.begin(), QMSegAtIt = QMSegAtoms.begin();
                    MDSegAtIt < MDSegAtoms.end();
                    MDSegAtIt++, QMSegAtIt++)
            {
                new_res_num = (*MDSegAtIt)->getResnr();
                if ( new_res_num != old_res_num)
                {
                    old_res_num = new_res_num;
                    res_name = (*MDSegAtIt)->getResname();
                    
                    QMfragment = _QMtop.AddFragment(res_name);
                    QMmolecule->AddFragment(QMfragment);
                    (*QMSegIt)->AddFragment(QMfragment);
                    
                    QMfragment->setTopology(&_QMtop);
                    QMfragment->setMolecule(QMmolecule);
                    QMfragment->setSegment(*QMSegIt);
                    
                    (*QMSegAtIt)->setFragment(QMfragment);
                    QMfragment->AddAtom(*QMSegAtIt);                    
                }
                else
                {
                    (*QMSegAtIt)->setFragment(QMfragment);
                    QMfragment->AddAtom(*QMSegAtIt);    
                }
            }
        }
    }
    else{
        error1("... ... Number of MD atoms is different from QM.\n"
               "... ... If it's the case of reduced molecule, "
                        " I need a map.\n"
               " ... ... Tags: map\n"
               " ... ... NOT IMPLEMENTED\n");
    }
}

void PDB2Map::readPDB(){
   cout << endl << "... ... Assuming: PDB for MD.";
    
    // set molecule >> segment >> fragment
    // reconnect them all
    Topology * _topPtr = 0;
    _topPtr = &_MDtop;
    
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
    int _newResNum = 0;
    bool warning_showed = false;

    while ( std::getline(_file, _line,'\n') ){
        if(     boost::find_first(_line, "ATOM"  )   || 
                boost::find_first(_line, "HETATM")      
                ){
            
            if ( (!_has_xyz && !warning_showed) || _line.size()<78 ){
               cout << endl << "... ... WARNING: No chemical elements in PDB!\n"
                            << "... ... Expect: empty slots \n"
                            << "in <qmatoms> and <multipoles>, "
                               "zeros in <weights>.\n"
                            << "... ... To add chemical symbols use: "
                               "editconf (GROMACS), babel, "
                               "(hands+pdb format)";                   
               warning_showed = true;
            }
            std::cout << "The size of str is " << _line.size() << " bytes.\n";
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
    
    // if all was read OK and
    // no XYZ file mentioned for QM top
    // make MD topology convertable to QM top
    if (!_has_xyz){_can_convert_md2qm = true;}
    
    return;
}

void PDB2Map::readGRO(){
    cout << endl << "... ... Assuming: GRO for MD.";

    // set molecule >> segment >> fragment
    // reconnect them all
    Topology * _topPtr = 0;
    _topPtr = &_MDtop;
    
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
            int _resNumInt(0);//,_atNumInt(0);
            double _xd(0),_yd(0),_zd(0);
            try
            {
                _resNumInt = boost::lexical_cast<int>(_resNum);
		//_atNumInt  = boost::lexical_cast<int>(_atNum);

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

void PDB2Map::readXYZ(){
    cout << endl << "... ... Assuming: XYZ for QM.";
    
    // set molecule >> segment >> fragment
    // reconnect them all
    Topology * _topPtr = 0;
    _topPtr = &_QMtop;
    
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
    
//    Fragment * _fragPtr = 0;
//    // direct
//    _fragPtr = _topPtr->AddFragment("F1");
//               _molPtr->AddFragment(_fragPtr);
//               _segPtr->AddFragment(_fragPtr);
//                // inverse
//              _fragPtr->setTopology(_topPtr);
//              _fragPtr->setMolecule(_molPtr);
//              _fragPtr->setSegment(_segPtr);
    
    // try: read xyz file
    std::ifstream _file( _input_xyz.c_str());
    if (_file.fail()) {
        error1( "... ... Can not open: " + _input_xyz + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    }
    else{
        cout << endl << 
                ("... ... File " + _input_xyz + ""
                 " was opened successfully.\n");
    }
    
    // read XYZ line by line
    string _line;
    
    // XYZ: first two lines are tech specs -> ignore them
    // XYZ check: if first line can cast to int, then ok
    try
    {   
        // first line, number of atoms in XYZ
        std::getline(_file, _line,'\n');
        ba::trim(_line);
        //int numXYZatoms = boost::lexical_cast<double>(_line);
    }
    catch(boost::bad_lexical_cast &)
    {
        error1( "... ... Bad XYZ file format!\n"
                "... ... First two line must contain technical specs.\n");
    }
    
    // ignore second line, it's a comment
    std::getline(_file, _line,'\n');
    
    while ( std::getline(_file, _line,'\n') ){

        // tokenize wrt space (free format)
        Tokenizer tokLine( _line, " ");
        vector<string> vecLine;
        tokLine.ToVector(vecLine);
        
        if (vecLine.size()!=4){
            error1("... ... Bad coord line in XYZ. Fix your XYZ file!\n");
        }
        
        string _atName     (vecLine[0]); // str,  Atom name
        string _x          (vecLine[1]); // 
        string _y          (vecLine[2]); // 
        string _z          (vecLine[3]); // 
        
        // try transform xyz coords to double
        double _xd(0),_yd(0),_zd(0);
        try{
            _xd = boost::lexical_cast<double>(_x);
            _yd = boost::lexical_cast<double>(_y);
            _zd = boost::lexical_cast<double>(_z);
        }
        catch(boost::bad_lexical_cast &)
        {
                 error1( "... ... Can't make numbers from strings.\n"
                         "... ... I don't like this string: \n\n"
                         "" + _line + "\n\n"
                         "... ... Check if coords are numbers!\n");
        }
        vec r(_xd , _yd , _zd);
        
        // set atom
        // reconnect to topology, molecule, segment, fragment
        Atom * _atmPtr = 0;
        // direct
        _atmPtr = _topPtr->AddAtom(_atName);
                _molPtr->AddAtom(_atmPtr);
                 _segPtr->AddAtom(_atmPtr);
//                _fragPtr->AddAtom(_atmPtr);
                    // inverse
                    _atmPtr->setTopology(_topPtr);
                    _atmPtr->setMolecule(_molPtr);        
                    _atmPtr->setSegment(_segPtr);
//                    _atmPtr->setFragment(_fragPtr);
        
        // set atom name, position
        _atmPtr->setElement(_atName);
        _atmPtr->setPos(r);
    }
    
    return;
}

void PDB2Map::topMdQm2xml(){
    cout << endl << "... ... (A)merging XML from MD and QM topologies.";
    if ( !_has_qm ) 
    {
        error1("... ... Error from topMdQm2xml(). QM top is missing.");    
    }
    else if (!_has_md)
    {
        error1("... ... Error from topMdQm2xml(). MD top is missing.");   
    }
    Molecule * MDmolecule = _MDtop.getMolecule(1);
    Molecule * QMmolecule = _QMtop.getMolecule(1);
    
    //
    // xml stuff
    //
    
    Property record;
    Property *ptopology_p = &record.add("topology","");
    Property *pmolecules_p = &ptopology_p->add("molecules","");
    Property *pmolecule_p = &pmolecules_p->add("molecule","");
    pmolecule_p->add("name","random_molecule_name");
    pmolecule_p->add("mdname","name_from_topol.top");
    Property *psegments_p = &pmolecule_p->add("segments","");
    Property *psegment_p = &psegments_p->add("segment","");
    psegment_p->add("name","random_segment_name");
    
    // qc data
    psegment_p->add("qmcoords","QC_FILES/your_file_with.xyz");
    psegment_p->add("orbitals","QC_FILES/your_file_with.fort7");
    psegment_p->add("basisset","INDO");
    
    // hole
    psegment_p->add("torbital_h","number_of_electrons,150");
    psegment_p->add("U_cC_nN_h","0.0000");
    psegment_p->add("U_nC_nN_h","0.0000");
    psegment_p->add("U_cN_cC_h","0.0000");
    
    // electron
    psegment_p->add("torbital_e","number_of_electrons+1,151");
    psegment_p->add("U_cC_nN_e","0.0000");
    psegment_p->add("U_nC_nN_e","0.0000");
    psegment_p->add("U_cN_cC_e","0.0000");

    // mps entry
    psegment_p->add("multipoles_n","MP_FILES/file_with.mps");
    psegment_p->add("multipoles_h","MP_FILES/file_with.mps");
    psegment_p->add("multipoles_e","MP_FILES/file_with.mps");
    
    psegment_p->add("map2md","0");
    
    // main body 
    Property *pfragments_p = &psegment_p->add("fragments","");
                                        
    vector < Segment * > allMdSegments = MDmolecule->Segments();
    vector < Segment * > allQmSegments = QMmolecule->Segments();
  
    vector < Segment * >::iterator segMdIt;
    vector < Segment * >::iterator segQmIt;
                                        
    for ( segMdIt = allMdSegments.begin(), 
                segQmIt = allQmSegments.begin();
            
          (allMdSegments.size() > allQmSegments.size()) ? 
              segMdIt < allMdSegments.end() :
              segQmIt < allQmSegments.end();
            
          segMdIt++, segQmIt++)
    {
        
        vector < Fragment * > allMdFragments = (*segMdIt)->Fragments();
        vector < Fragment * > allQmFragments = (*segQmIt)->Fragments();
                                        
        vector < Fragment * >::iterator fragMdIt;
        vector < Fragment * >::iterator fragQmIt;
                                        
        for ( fragMdIt = allMdFragments.begin() , 
                fragQmIt = allQmFragments.begin();
                
              (allMdFragments.size() > allQmFragments.size()) ?
                  fragMdIt < allMdFragments.end() :
                  fragQmIt < allQmFragments.end();
                
              fragMdIt++,fragQmIt++ )
        {
            string mapName;            
            stringstream mapMdAtoms;
            stringstream mapQmAtoms;
            stringstream mapMpoles;
            stringstream mapWeight;
            stringstream mapFrame;

            mapName      = (*fragMdIt)->getName() ;
            
            unsigned int localCounter = 0;
            vector < Atom * > allMdAtoms = (*fragMdIt)->Atoms();
            vector < Atom * > allQmAtoms = (*fragQmIt)->Atoms();

            vector < Atom * >::iterator atomMdIt;
            vector < Atom * >::iterator atomQmIt;
            for ( atomMdIt = allMdAtoms.begin(),
                    atomQmIt = allQmAtoms.begin();
                    
                  (allMdAtoms.size() > allQmAtoms.size()) ? 
                      atomMdIt < allMdAtoms.end() :
                      atomQmIt < allQmAtoms.end();
                  
                  atomMdIt++, atomQmIt++ )
            {
                if (atomMdIt < allMdAtoms.end())
                {
                        mapMdAtoms << boost::format("%=13s") % 
                                (boost::format("%s:%s:%s") 
                                   % (*atomMdIt)->getResnr()
                                   % (*atomMdIt)->getResname()
                                   % (*atomMdIt)->getName()).str() ;
                        
                                if (el2mass.find((*atomQmIt)->getElement()) 
                                        != el2mass.end())
                                {
                        mapWeight << boost::format("%=13i")
                                           % el2mass[(*atomQmIt)->getElement()];
                                }
                                else
                                {
                        mapWeight << boost::format("%=13i") % " ";        
                                }
                }
                else
                {
                        mapMdAtoms << boost::format("%=13s") %
                                (boost::format("%s:%s:%s") 
                                        % " " % " " % " " ).str();
                        
                        mapWeight << boost::format("%=13i") % " " ;
                }
                
                if (atomQmIt < allQmAtoms.end())
                {
                        mapQmAtoms << boost::format("%=13s") %
                                (boost::format("%1%:%2% ") 
                                   % (*atomQmIt)->getId()
                                   % (*atomQmIt)->getElement()).str() ;
                }
                else
                {
                        mapQmAtoms << boost::format("%=13s") %
                                (boost::format("%1%:%2% ") 
                                        % " " % " " ).str();
                }
                
                if (localCounter < 3 && localCounter < allMdAtoms.size())
                {
                        mapFrame << boost::format("%=5i")
                                   % (*atomMdIt)->getId();
                }
                localCounter++;
            }
            
            mapMpoles << " " << mapQmAtoms.str();
            
            Property *pfragment_p  = &pfragments_p->add("fragment","");
            pfragment_p->add("name", mapName);
            pfragment_p->add("mdatoms", mapMdAtoms.str());
            pfragment_p->add("qmatoms", mapQmAtoms.str());
            pfragment_p->add("mpoles",  mapMpoles.str());
            pfragment_p->add("weights", mapWeight.str());
            pfragment_p->add("localframe", mapFrame.str());

         }
    }
    
//    cout << endl << setlevel(1) << XML << record;
        
    ofstream outfile( _output_xml.c_str() );
    
    PropertyIOManipulator XML(PropertyIOManipulator::XML,1,""); 
    outfile << XML << record;
    outfile.close();
    
    
    return;
}

}}


#endif
