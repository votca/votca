#ifndef _PDB2Map_H
#define _PDB2Map_H


#include <votca/ctp/topology.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/logger.h>
#include <boost/algorithm/string.hpp>
#include <votca/tools/vec.h>
#include <boost/format.hpp>


namespace votca { namespace ctp {
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
    void compatibilityQM2MD();
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
    
    // good guys
    _has_pdb = false;
    _has_gro = false;
    _has_xyz = false;
    
    _has_qm  = false;
    _has_md  = false;
    
    _can_convert_md2qm = false;
    
    // read options    
    string key = "options.pdb2map.";
    
    // find PDB, then GRO, then error
    if ( options->exists(key+"pdb") ){
        _input_pdb      = options->get(key+"pdb").as<string> ();
        _has_pdb        = true;
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
    
    // find XYZ 
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
    compatibilityQM2MD();
    topMdQm2xml();
    
//    LOG( logINFO, _log ) << "Reading from: " << _input_file << flush;    
//    cout << _log;
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
        error1("No good MD topology. I stop.");
    }
    
    if (_has_xyz){
        readXYZ();
        _has_qm = true;
    }
    else if (_can_convert_md2qm){
        _has_qm = true;
    }
    else{
        error1("No good QM topology. I stop.");
    }
}

void PDB2Map::compatibilityQM2MD(){
    int numMDatoms = _MDtop.getMolecule(1)->NumberOfAtoms();
    int numQMatoms = _QMtop.getMolecule(1)->NumberOfAtoms();
    
    _QM2MDcompatible = (numMDatoms == numQMatoms) ? true : false;
    
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
            Fragment * QMfragment;
            
            vector<Atom*> MDSegAtoms = (*MDSegIt)->Atoms();
            vector<Atom*> QMSegAtoms = (*QMSegIt)->Atoms();
            
            vector<Atom*>::iterator MDSegAtIt;
            vector<Atom*>::iterator QMSegAtIt;
            
            int old_res_num = -1;
            int new_res_num;
            string res_name = "bad_wolf";
            
            for(MDSegAtIt = MDSegAtoms.begin(), QMSegAtIt = QMSegAtoms.begin();
                    MDSegAtIt < MDSegAtoms.end();
                    MDSegAtIt++, QMSegAtIt++)
            {
                new_res_num = (*MDSegAtIt)->getResnr();
                if (new_res_num != old_res_num)
                {
                    old_res_num = new_res_num;
                    res_name = (*MDSegAtIt)->getResname();
                    
                    QMfragment = _QMtop.AddFragment(res_name);
                    QMfragment->setTopology(&_QMtop);
                    
                    QMmolecule->AddFragment(QMfragment);
                    QMfragment->setMolecule(QMmolecule);
                    
                    (*QMSegIt)->AddFragment(QMfragment);
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
        error1("\n... ... Number of MD atoms is different from QM."
               "\n... ... If it's the case of reduced molecule, "
                        " I need a map."
               "\n ... ... Tags: map");
    }
}

void PDB2Map::readPDB(){
    
    cout << endl << "... ... Assuming: PDB for MD. Read.";
    
    // make molecule and segment in molecule
    Molecule * newMolecule = _MDtop.AddMolecule("newMolecule");
    newMolecule->setTopology(&_MDtop);
    
    Segment  * newSegment  = _MDtop.AddSegment ("newSegment");
    newSegment->setTopology(&_MDtop);
    
    newMolecule->AddSegment(newSegment);
    newSegment->setMolecule(newMolecule);

    // reading from PDB file and creating topology
    std::ifstream _file( _input_pdb.c_str() );
    
    if (!_file.is_open()) {
       error1(  "... ... Bad file: " + _input_pdb + \
              "\n... ... Does it exist? Bad name?");
    }
    
    string _line;

    int _atom_id = 0;
    int _newResNum = 0;
    bool chem_message_showed = false;

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
            
            if (_atElement.empty() && !chem_message_showed && !_has_xyz ){
                cout << endl << "... ... *** No chemical elements in PDB!"
                        << endl << "... ... *** Expect: empty slots "
                        << "in <qmatoms> and <multipoles>, zeros in <weights>.";
                chem_message_showed = true;
            }
            else
            {
              _can_convert_md2qm = true;
            }
            

            double _xd = boost::lexical_cast<double>(_x);
            double _yd = boost::lexical_cast<double>(_y);
            double _zd = boost::lexical_cast<double>(_z);
            int _resNumInt = boost::lexical_cast<int>(_resNum);

            vec r(_xd , _yd , _zd);

            Atom * newAtom = _MDtop.AddAtom(_atName);
            newAtom->setTopology(&_MDtop);

            newAtom->setResnr        (_resNumInt);
            newAtom->setResname      (_resName);
            newAtom->setElement      (_atElement);
            newAtom->setPos          (r);

            newMolecule->AddAtom(newAtom);
            newAtom->setMolecule(newMolecule);        
            
            newSegment->AddAtom(newAtom);
            newAtom->setSegment(newSegment);

            Fragment * newFragment;
            if ( _newResNum != _resNumInt ){

                _newResNum = _resNumInt;
                string _newResName = _resName+'_'+_resNum;
                
                newFragment = _MDtop.AddFragment(_newResName);
                newFragment->setTopology(&_MDtop);

                newMolecule->AddFragment(newFragment);
                newFragment->setMolecule(newMolecule);

                newSegment->AddFragment(newFragment);
                newFragment->setSegment(newSegment);
            }

            newFragment->AddAtom(newAtom);
            newAtom->setFragment(newFragment);
        }
    }
    
    return;
}

void PDB2Map::readGRO(){
    
    cout << endl << "... ... Assuming: GRO for MD. Read.";

    // make molecule and segment in molecule
    Molecule * newMolecule = _MDtop.AddMolecule("newMolecule");
    newMolecule->setTopology(&_MDtop);
    
    Segment  * newSegment  = _MDtop.AddSegment ("newSegment");
    newSegment->setTopology(&_MDtop);
    
    newMolecule->AddSegment(newSegment);
    newSegment->setMolecule(newMolecule);

    // reading from GRO file and creating topology
    std::ifstream _file( _input_gro.c_str() );

    if (!_file.is_open()) {
        cout << endl;
        throw runtime_error("Bad file: " + _input_pdb + \
                            "\n... ... Does it exist? Bad name?");
    }
    
    string _line;

    int _atom_id = 0;
    int _newResNum = 0;
    
    // ignore first 2 lines - as in GRO format
    std::getline(_file, _line,'\n');
    std::getline(_file, _line,'\n');
    
    ba::trim(_line);
    int atom_num = boost::lexical_cast<int>(_line);
    int counter = 0;
    
    while ( std::getline(_file, _line,'\n') ){
        if (counter < atom_num){
            
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
            
            double _xd = boost::lexical_cast<double>(_x);
            double _yd = boost::lexical_cast<double>(_y);
            double _zd = boost::lexical_cast<double>(_z);
            
            vec r(_xd , _yd , _zd);
            
            int _resNumInt = boost::lexical_cast<int>(_resNum);
            
            Atom * newAtom = _MDtop.AddAtom(_atName);
            newAtom->setTopology(&_MDtop);

            newAtom->setResnr        (_resNumInt);
            newAtom->setResname      (_resName);
            newAtom->setPos          (r);

            newMolecule->AddAtom(newAtom);
            newAtom->setMolecule(newMolecule);        
            
            newSegment->AddAtom(newAtom);
            newAtom->setSegment(newSegment);

            Fragment * newFragment;
            if ( _newResNum != _resNumInt ){

                _newResNum = _resNumInt;
                string _newResName = _resName+'_'+_resNum;
                
                newFragment = _MDtop.AddFragment(_newResName);
                newFragment->setTopology(&_MDtop);

                newMolecule->AddFragment(newFragment);
                newFragment->setMolecule(newMolecule);

                newSegment->AddFragment(newFragment);
                newFragment->setSegment(newSegment);
            }

            newFragment->AddAtom(newAtom);
            newAtom->setFragment(newFragment);
        }
        counter++;
    }
    
    return;
}

void PDB2Map::readXYZ(){
    cout << endl << "... ... Assuming: XYZ for QM. Read.";
    
    // make molecule and segment
    Molecule * newMolecule = _QMtop.AddMolecule("newMolecule");
    newMolecule->setTopology(&_QMtop);
    
    Segment  * newSegment  = _QMtop.AddSegment ("newSegment");
    newSegment->setTopology(&_QMtop);
    newSegment->setMolecule(newMolecule);
    newMolecule->AddSegment(newSegment);
    
//    Fragment * newFragment = _QMtop.AddFragment("newFragment");
//    newFragment->setTopology(&_QMtop);
//    newMolecule->AddFragment(newFragment);
//    newFragment->setMolecule(newMolecule);
//    newSegment->AddFragment(newFragment);
//    newFragment->setSegment(newSegment);
    
    // reading from PDB file and creating topology
    std::ifstream _file( _input_xyz.c_str() );
    
    if (!_file.is_open()) {
        error1(   "... ... Bad file: " + _input_pdb + \
                "\n... ... Does it exist? Bad name?");
    }    
       
    string _line;
    
    // ignoring first 2 lines - as in XYZ format
    std::getline(_file, _line,'\n');
    std::getline(_file, _line,'\n');
    
    while ( std::getline(_file, _line,'\n') ){
 
        Tokenizer tokLine( _line, " ");
        vector<string> vecLine;
        tokLine.ToVector(vecLine);
    
        string _atName     (vecLine[0]); // str,  Atom name
        string _x          (vecLine[1]); // float 8.3 ,x
        string _y          (vecLine[2]); // float 8.3 ,y
        string _z          (vecLine[3]); // float 8.3 ,z

        double _xd = boost::lexical_cast<double>(_x);
        double _yd = boost::lexical_cast<double>(_y);
        double _zd = boost::lexical_cast<double>(_z);

        vec r(_xd , _yd , _zd);
            
        Atom * newAtom = _QMtop.AddAtom(_atName);
        newAtom->setElement(_atName);
        newAtom->setTopology(&_QMtop);
        newAtom->setPos(r);

        newMolecule->AddAtom(newAtom);
        newAtom->setMolecule(newMolecule);        

        newSegment->AddAtom(newAtom);
        newAtom->setSegment(newSegment);
        
//        newFragment->AddAtom(newAtom);
//        newAtom->setFragment(newFragment);
    }
    
    return;
}

void PDB2Map::topMdQm2xml(){
    cout << endl << "... ... (A)merging XML from MD and QM topologies.";
    
    Molecule * MDmolecule = _MDtop.getMolecule(1);
    Molecule * QMmolecule = _QMtop.getMolecule(1);
    // xml stuff
    
    Property record;
    Property *ptopology_p = &record.add("topology","");
    Property *pmolecules_p = &ptopology_p->add("molecules","");
    Property *pmolecule_p = &pmolecules_p->add("molecule","");
    pmolecule_p->add("name","MOLECULE_NAME");
    pmolecule_p->add("mdname","Other");
    Property *psegments_p = &pmolecule_p->add("segments","");
    Property *psegment_p = &psegments_p->add("segment","");
    psegment_p->add("name","SEGMENT_NAME");
    psegment_p->add("qmcoords","QC_FILES/your_file_with.xyz");
    psegment_p->add("orbitals","QC_FILES/your_file_with.fort7");
    psegment_p->add("basisset","INDO");
    psegment_p->add("torbital_h","NUMBER");
    psegment_p->add("U_cC_nN_h","NUMBER");
    psegment_p->add("U_nC_nN_h","NUMBER");
    psegment_p->add("U_cN_cC_h","NUMBER");
    psegment_p->add("multipoles_n","MP_FILES/your_file_with.mps");
    psegment_p->add("multipoles_h","your_file_with.mps");
    psegment_p->add("map2md","1");
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
            
            int localCounter = 0;
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
    outfile << setlevel(1) << XML << record;
    outfile.close();

    return;
}

}}


#endif