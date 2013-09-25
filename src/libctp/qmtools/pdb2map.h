#ifndef _PDB2Map_H
#define _PDB2Map_H


#include <votca/ctp/topology.h>
#include <votca/ctp/atom.h>
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

    void   Initialize(Property *options);
    bool   Evaluate();
    
    void readPDB();
    void readGRO();
    void readXYZ();
    void topMdQm2xml();
 

private:
    string      _input_pdb;
    string      _input_gro;
    string      _input_xyz;
    string      _output_xml;
    
    bool        _has_pdb;
    bool        _has_gro;
    bool        _has_xyz;    
    
    string      _print_pdb;
    
    string      _input_file;
    string      _output_file;
    Topology    _MDtop;
    Topology    _QMtop;

   // element:mass map
    map <string,int> el2mass;

};

void PDB2Map::Initialize(Property* options) 
{   

    // fill out periodic table
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

    _has_pdb = false;
    _has_gro = false;
    _has_xyz = false;
    
    // options reading    
    string key = "options.pdb2map.";
    
    // old options
    _input_file  = options->get(key+"file").as<string> ();
    _output_file = options->get(key+"outfile").as<string> ();
    
    // new options
    // find PDB or GRO
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
        throw runtime_error("... ... No MD(PDB,GRO) file provided. Tags: pdb, gro");
    }
    
    // find XYZ 
    if ( options->exists(key+"xyz") ){
        _input_xyz      = options->get(key+"xyz").as<string> ();
        _has_xyz        = true;
        cout << endl << "... ... XYZ input specified: \t" << _input_xyz;
    }
    else{
        throw runtime_error("... ... No QM(XYZ) file provided. Tags: xyz");
    }
    
    // find XML or generate it
    if ( options->exists(key+"xml") ){
        _output_xml = options->get(key+"xml").as<string> ();
        cout << endl << "... ... XML output specified: \t" << _output_xml;
        }
        else{
                _output_xml = "system_output.xml";
                cout << endl << "... ... No XML output specified.";
                cout << endl << "... ... Default XML is: \t" << _output_xml;
        }
}



bool PDB2Map::Evaluate() {
    
    if ( _has_pdb ){
        readPDB();
    }
    else if ( _has_gro ){
        readGRO();
    }

    if ( _has_xyz ){
        readXYZ();
    }
    
    topMdQm2xml();
    
//    LOG( logINFO, _log ) << "Reading from: " << _input_file << flush;    
//    cout << _log;
    
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
    string _line;

    int _atom_id = 0;
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

    // reading from PDB file and creating topology
    std::ifstream _file( _input_gro.c_str() );
    string _line;

    int _atom_id = 0;
    int _newResNum = 0;
    
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
//
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
    
    // make molecule and segment in molecule
    Molecule * newMolecule = _QMtop.AddMolecule("newMolecule");
    newMolecule->setTopology(&_QMtop);

    Segment  * newSegment  = _QMtop.AddSegment ("newSegment");
    newSegment->setTopology(&_QMtop);
    newSegment->setMolecule(newMolecule);
    newMolecule->AddSegment(newSegment);
    
    Fragment * newFragment = _QMtop.AddFragment("newFragment");
    newFragment->setTopology(&_QMtop);
    newMolecule->AddFragment(newFragment);
    newFragment->setMolecule(newMolecule);
    newSegment->AddFragment(newFragment);
    newFragment->setSegment(newSegment);

    // reading from PDB file and creating topology
    std::ifstream _file( _input_xyz.c_str() );
    string _line;
    
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
        newAtom->setTopology(&_QMtop);

        newAtom->setPos          (r);

        newMolecule->AddAtom(newAtom);
        newAtom->setMolecule(newMolecule);        

        newSegment->AddAtom(newAtom);
        newAtom->setSegment(newSegment);
        
        newFragment->AddAtom(newAtom);
        newAtom->setFragment(newFragment);
    }    
    
    return;
}

void PDB2Map::topMdQm2xml(){
    cout << endl << "... ... Amerging XML from MD and QM topologies.";
    
    Molecule * newMolecule = _MDtop.getMolecule(1);
    
    // iterating over thins, making system.xml
    stringstream map2pdbFile;
    
    Property record;
    Property *pfragment = &record.add("topology","");
    pfragment->add("molecules.new", "");
    pfragment->add("molecules", "");
    pfragment->add("molecule", "");
    pfragment->add("name", "MOLECULE_NAME");
    pfragment->add("molecule", "MOLECULE_NAME");
    pfragment->add("molecule", "");
    pfragment->add("molecule", "");
    pfragment->add("molecule", "");
    pfragment->add("molecule", "");
    
    string headlines = 
                  "<topology>\n"
                "\t<molecules>\n"
                "\t<molecule>\n"
              "\t\t<name>MOLECULE_NAME</name>\n"
              "\t\t<mdname>Other</mdname>\n"
              "\t\t<segments>\n"
              "\t\t<segment>\n"
            "\t\t\t<name>SEGMENT_NAME</name>\n"
            "\t\t\t<qmcoords>QC_FILES/your_file_with.xyz</qmcoords>\n\n"
            "\t\t\t<orbitals>QC_FILES/your_file_with.fort7</orbitals>\n"
            "\t\t\t<basisset>INDO</basisset>\n"
            "\t\t\t<torbital_h>NUMBER</torbital_h>\n\n"
            "\t\t\t<U_cC_nN_h>NUMBER</U_cC_nN_h>\n"
            "\t\t\t<U_nC_nN_h>NUMBER</U_nC_nN_h>\n"
            "\t\t\t<U_cN_cC_h>NUMBER</U_cN_cC_h>\n\n"
            "\t\t\t<multipoles_n>MP_FILES/your_file_with.mps</multipoles_n>\n"
            "\t\t\t<multipoles_h>MP_FILES/your_file_with.mps</multipoles_h>\n"
            "\t\t\t<map2md>1</map2md>\n\n"
            "\t\t\t<fragments>\n";
    
    map2pdbFile << headlines;
    
    
    vector < Segment * > allSegments = newMolecule->Segments();
    vector < Segment * >::iterator segIt;
    for (       segIt = allSegments.begin();
                segIt < allSegments.end();
                segIt++){
        
        vector < Fragment * > allFragments = (*segIt)->Fragments();
        vector < Fragment * >::iterator fragIt;
        for (   fragIt = allFragments.begin();
                fragIt < allFragments.end();
                fragIt++){
            
            string mapName      = (*fragIt)->getName() ;
            stringstream mapMdAtoms;
            stringstream mapQmAtoms;
            stringstream mapWeight;
            stringstream mapFrame;
            
            int localCounter = 0;
            vector < Atom * > allAtoms = (*fragIt)->Atoms();
            vector < Atom * >::iterator atomIt;
            for (       atomIt = allAtoms.begin();
                        atomIt < allAtoms.end();
                        atomIt++){
                
               mapMdAtoms << boost::format("%=13s") % 
                                (boost::format("%s:%s:%s") 
                                   % (*atomIt)->getResnr()
                                   % (*atomIt)->getResname()
                                   % (*atomIt)->getName()).str();
                mapQmAtoms << boost::format("%=13s") %
                                (boost::format("%1%:%2% ") 
                                   % (*atomIt)->getId()
                                   % (*atomIt)->getElement()).str();
                mapWeight << boost::format("%=13i") 
                                   % el2mass[(*atomIt)->getElement()];
                
                if (localCounter < 3){
                        mapFrame << boost::format("%=5i")
                                   % (*atomIt)->getId();
                }
                localCounter++;
            }
            

//            pfragment->add("name", mapName);
//            pfragment->add("mdatoms", mapMdAtoms.str());
//            pfragment->add("qmatoms", mapQmAtoms.str());
//            pfragment->add("mpoles", mapQmAtoms.str());
//            pfragment->add("weights", mapWeight.str());
//            pfragment->add("localframe", mapFrame.str());
                    

            
            map2pdbFile 
                 << "\t\t\t\t<fragment>\n"
                 << "\t\t\t\t\t<name> "     << mapName          << " </name>\n"
                 << "\t\t\t\t\t<mdatoms>"   << mapMdAtoms.str() << "</mdatoms>\n"
                 << "\t\t\t\t\t<qmatoms>"   << mapQmAtoms.str() << "</qmatoms>\n"
                 << "\t\t\t\t\t<mpoles> "   << mapQmAtoms.str() << "</mpoles>\n"
                 << "\t\t\t\t\t<weights>"   << mapWeight.str()  << "</weights>\n"
                 << "\t\t\t\t\t<localframe>"<< mapFrame.str()   << "</localframe>\n"   
                 << "\t\t\t\t</fragment>\n";
         }
    }
    
    string tailLine =
    "\t\t\t</fragments>\n"
    "\t\t\t</segment>\n"              
      "\t\t</segments>\n"
        "\t</molecule>\n"    
        "\t</molecules>\n"
          "</topology>\n";    
    map2pdbFile << tailLine;
    
    cout << endl << setlevel(1) << XML << record;
        
    ofstream outfile( _output_file.c_str() );
    outfile << map2pdbFile.str();
    outfile.close();

    return;
}

}}


#endif