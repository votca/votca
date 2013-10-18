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

 

private:
    string      _input_file;
    string      _output_file;
    Topology    _top;
    Logger      _log;

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

    
    // options reading
    _input_file  = options->get("options.pdb2map.file").as<string> ();
    _output_file = options->get("options.pdb2map.outfile").as<string> ();
    
    //cout << "\n... ... " <<"Reading from: " << _input_file << "\n\n";

    
}



bool PDB2Map::Evaluate() {
    
    LOG( logINFO, _log ) << "Reading from: " << _input_file << flush;
    
    // read strings to vector
    
    std::ifstream _file( _input_file.c_str() );
    
    vector<string>      strVec;
    string _line;
    // push lines to vector + only ATOM, HETATM
    while ( std::getline(_file, _line,'\n') ){
        if(boost::find_first(_line, "ATOM") || 
                boost::find_first(_line, "HETATM")){
            strVec.push_back(_line);
        }
    }
    
    // topology insanity
    
    Molecule * newMolecule = _top.AddMolecule("newMolecule");
    Segment  * newSegment  = _top.AddSegment ("newSegment");
    
    newMolecule->AddSegment(newSegment);
    newSegment->setMolecule(newMolecule);
            
 
    // iterate over lines, create atoms
    
    int _atom_id = 0;
    int _newResNum = 0;
    for (vector<string>::iterator p = strVec.begin() ; 
                        p != strVec.end(); ++p)
    {
        
        //      according to PDB format
        string _recType    (*p,( 1-1),6); // str,  "ATOM", "HETATM"
        string _atNum      (*p,( 7-1),6); // int,  Atom serial number
        string _atName     (*p,(13-1),4); // str,  Atom name
        string _atAltLoc   (*p,(17-1),1); // char, Alternate location indicator
        string _resName    (*p,(18-1),4); // str,  Residue name
        string _chainID    (*p,(22-1),1); // char, Chain identifier
        string _resNum     (*p,(23-1),4); // int,  Residue sequence number
        string _atICode    (*p,(27-1),1); // char, Code for insertion of res
        string _x          (*p,(31-1),8); // float 8.3 ,x
        string _y          (*p,(39-1),8); // float 8.3 ,y
        string _z          (*p,(47-1),8); // float 8.3 ,z
        string _atOccup    (*p,(55-1),6); // float  6.2, Occupancy
        string _atTFactor  (*p,(61-1),6); // float  6.2, Temperature factor
        string _segID      (*p,(73-1),4); // str, Segment identifier
        string _atElement  (*p,(77-1),2); // str, Element symbol
        string _atCharge   (*p,(79-1),2); // str, Charge on the atom
        
        
        //              I print out formated pdb if print_us = 1
        int print_us = 0;
        if (print_us == 1)
        {
            string _delim = "|";

            cout << _recType   << _delim;
            cout << _atNum     << _delim;
            cout << _atName    << _delim;
            cout << _atAltLoc  << _delim;
            cout << _resName   << _delim;
            cout << _chainID   << _delim;
            cout << _resNum    << _delim;
            cout << _atICode   << _delim;
            cout << _x         << _delim;
            cout << _y         << _delim;
            cout << _z         << _delim;
            cout << _atOccup   << _delim;
            cout << _atTFactor << _delim;
            cout << _segID     << _delim;
            cout << _atElement << _delim;
            cout << _atCharge  <<   endl; 
        }
        
        
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
        
        Atom *newAtom = _top.AddAtom(_atName);
        
        newAtom->setResnr        (_resNumInt);
        newAtom->setResname      (_resName);
        newAtom->setElement      (_atElement);
        newAtom->setPos          (r);
        
        newAtom->setMolecule(newMolecule);        
        newAtom->setSegment(newSegment);
        
        newMolecule->AddAtom(newAtom);
        newSegment->AddAtom(newAtom);
        
        Fragment * newFragment;
        if ( _newResNum != _resNumInt ){
            _newResNum = _resNumInt;
            string _newResName = _resName+'_'+_resNum;
            newFragment = _top.AddFragment(_newResName);

            newFragment->setMolecule(newMolecule);
            newFragment->setSegment(newSegment);
            
            newMolecule->AddFragment(newFragment);
            newSegment->AddFragment(newFragment);
        }
        
        newFragment->AddAtom(newAtom);
    }
    
    // iterating over thins, making system.xml
    
    stringstream map2pdbFile;
    
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
            
            Property record;
         
            Property *pfragment = &record.add("fragment","");
            pfragment->add("name", mapName);
            pfragment->add("mdatoms", mapMdAtoms.str());
            pfragment->add("qmatoms", mapQmAtoms.str());
            pfragment->add("mpoles", mapQmAtoms.str());
            pfragment->add("weights", mapWeight.str());
            pfragment->add("localframe", mapFrame.str());

            votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");

            cout << iomXML << record;
            
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
    
    //cout << "\n... ... " <<"Writing into: " << _output_file << "\n\n";
    
    LOG( logINFO, _log ) << "Writing into: " << _output_file << flush;
        
    ofstream outfile( _output_file.c_str() );
    outfile << map2pdbFile.str();
    outfile.close();
    
    cout << _log;
    
    return true;
}


}}


#endif